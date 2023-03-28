//
// Created by Robert John Anderson on 2020-04-03.
//

#include "M7_lib/basis/Suites.h"
#include "Wavefunction.h"
#include "FciInitializer.h"
#include "HfExcitHists.h"


v_t<TableBase::Loc> wf::Vectors::setup() {
    v_t<TableBase::Loc> ref_locs;
    /*
     * create the reference MBF and add it to the walker table
     */
    buffered::Mbf ref_mbf(m_sector);
    mbf::set(ref_mbf, m_sector.particles(), m_opts.m_reference.m_mbf_init, 0ul);


    const auto permanitiators = m_opts.m_wavefunction.m_ci_permanitiator.m_enabled;

    if (permanitiators) {
        /*
         * using the "permanitiator" adaptation with low-rank CI information as a source
         */
        hf_excit_hist::initialize(*this, ref_mbf, m_opts.m_wavefunction.m_ci_permanitiator);
    }

    /*
     * insert reference MBF into the store table
     */
    const auto ref_loc = create_row_setup(0, ref_mbf);
    if (ref_loc.is_mine()) {
        auto ref_walker = m_store.m_row;
        ref_walker.jump(ref_loc.m_irec);
        for (uint_t ipart = 0ul; ipart < npart(); ++ipart) {
            set_weight(ref_walker, ipart, wf_t(m_opts.m_wavefunction.m_nw_init));
        }
        if (permanitiators) ref_walker.m_permanitiator.set();
    }

    for (auto ipart=0ul; ipart<npart(); ++ipart) ref_locs.push_back(ref_loc);

    if (m_opts.m_wavefunction.m_load.m_enabled) {
        // the wavefunction is to be loaded from HDF5 archive
    }
    else if (m_opts.m_wavefunction.m_fci_init) {
        // the wavefunction is to be initialized using exact eigenvectors from the Arnoldi method
        logging::info("Performing exact FCI initialization of wavefunctions");
        FciInitOptions fci_init_opts;
        fci_init_opts.m_nroot = this->nroot();
        fci_init(fci_init_opts);
    }
    return ref_locs;
}

wf::Vectors::Vectors(const conf::Document& opts, const Hamiltonian& ham):
    communicator::BasicSend<Walker, Spawn>(
        "wavefunction",
        // walker row:
        {
            ham.m_basis,
            opts.m_wavefunction.m_nroot,
            opts.m_av_ests.any_bilinears() ? 2ul:1ul, need_av_weights(opts)
        },
        opts.m_wavefunction.m_distribution,
        // store sizing
        {
            uint_t(opts.m_propagator.m_nw_target),
            opts.m_wavefunction.m_buffers.m_store_exp_fac
        },
        // send/recv row
        {ham.m_basis, need_send_parents(opts)},
        // send/recv sizing
        {
            std::max(10ul, uint_t(opts.m_propagator.m_nw_target * opts.m_propagator.m_tau_init)),
            opts.m_wavefunction.m_buffers.m_comm_exp_fac
        }
    ),
    m_opts(opts),
    m_ham(ham),
    m_sector(m_ham.m_basis, m_ham.default_particles(m_opts.m_particles)),
    m_format(m_store.m_row.m_weight.m_format),
    m_stats(m_format),
    m_refs(opts.m_reference, *this, setup()) {

    REQUIRE_TRUE(m_send_recv.recv().m_row.m_dst_mbf.belongs_to_row(), "row-field reference error");

    logging::info("Distributing wavefunction rows in {} block{}", m_dist.nblock(),
                  string::plural(m_dist.nblock()));
    refresh_all_hdiags();
    refresh_all_ref_conns();
}

void wf::Vectors::log_top_weighted(uint_t ipart, uint_t nrow) {
    buffered::Table<Walker> xr_gathered("global top weighted", m_store.m_row);
    {
        auto row1 = m_store.m_row;
        auto row2 = row1;
        weights_gxr_t gxr(row1.m_weight, row2.m_weight, true, true, ipart);
        gxr.find(nrow);
        gxr.gatherv(xr_gathered);
    }

    if (!mpi::i_am_root()) return;
    /*
     * the gathered rows (walkers) are globally maximal in occupation for component ipart, but they are simply laid
     * together by the gathering operation, and are not sorted internally. Here, that sorting operation is done
     */
    auto row1 = xr_gathered.m_row;
    auto row2 = row1;
    auto cmp_fn = [&](uint_t irow1, uint_t irow2){
        row1.jump(irow1);
        row2.jump(irow2);
        return std::abs(row1.m_weight[ipart]) > std::abs(row2.m_weight[ipart]);
    };

    quicksort::Sorter qs(cmp_fn);
    qs.reorder_sort(xr_gathered);

    auto& row = xr_gathered.m_row;
    v_t<strv_t> rows;
    rows.push_back({"", "many-body basis function", "walkers", "coefficient", "initiator", "energy", "semistoch", "MPI rank"});
    for (row.restart(); row; ++row) {
        rows.push_back({
            std::to_string(row.index()),
            row.m_mbf.to_string(),
            convert::to_string(row.m_weight[ipart], {true, 6}),
            convert::to_string(row.m_weight[ipart] / std::sqrt(m_stats.m_l2_norm_square.total()[ipart]), {false, 4}),
            convert::to_string(row.exceeds_initiator_thresh(ipart, m_opts.m_propagator.m_nadd) || row.m_permanitiator.get(0)),
            convert::to_string(row.m_hdiag[iroot_part(ipart)]),
            convert::to_string(bool(row.m_deterministic[iroot_part(ipart)])),
            convert::to_string(m_dist.irank(row.m_mbf))
        });
    }
    logging::info_table("Top-weighted WF elements for part "+std::to_string(ipart), rows, true, false, 1ul);
}

wf::Vectors::~Vectors() {
    for (uint_t ipart=0ul; ipart<npart(); ++ipart) log_top_weighted(ipart);
    if (m_opts.m_wavefunction.m_save.m_enabled) save();
}

void wf::Vectors::preserve_ref_weights(wf_comp_t mag) {
    if (m_ref_weights_preserved)
        return; // reference weights already being preserved
    for (uint_t ipart = 0ul; ipart < npart(); ++ipart) {
        if (!m_refs[ipart].is_mine()) continue;
        auto irec = m_refs[ipart].irec();
        auto row = m_store.m_row;
        row.jump(irec);
        const wf_t weight = row.m_weight[ipart];
        const auto fixed_weight = math::phase(weight) * mag;
        set_weight(row, ipart, fixed_weight);
    }
    m_ref_weights_preserved = true;
}

bool wf::Vectors::ref_weights_preserved() const {
    return m_ref_weights_preserved;
}

#if 0
void wf::Vectors::h5_write(const hdf5::NodeWriter& parent, str_t name) {
    auto field_names = []() -> strv_t {
        if (c_enable_fermions != c_enable_bosons) return {"mbf", "weight"};
        else if (c_enable_fermions != c_enable_bosons) return {"mbf (fermion)", "mbf (boson)", "weight"};
        return {};
    };
    m_store.save(parent, name, field_names(), true);
}

void wf::Vectors::h5_read(const hdf5::NodeReader& /*parent*/, const field::Mbf& /*ref*/,
                          str_t /*name*/) {
    m_store.clear();
    buffered::Table<Walker> m_buffer("", {m_store.m_row});
    m_buffer.push_back();
    RowHdf5Reader<Walker> row_reader(m_buffer.m_row, parent, name, h5_field_names());
    suite::Conns conn(m_sector.size());

    row_reader.restart();
    DEBUG_ASSERT_EQ(row_reader.m_weight.nelement(), m_format.m_nelement, "row reader has incompatible dimensionality");
    for (uint_t iitem = 0ul; iitem < row_reader.m_nitem; ++iitem) {
        row_reader.read(iitem);
        conn[ref].connect(ref, row_reader.m_mbf);
        bool ref_conn = ham::is_significant(ham.get_element(ref, conn[ref]));
        auto& walker = create_row(0ul, row_reader.m_mbf, ham.get_energy(row_reader.m_mbf), v_t<bool>(npart(), ref_conn));
        set_weight(walker, row_reader.m_weight);
    }
}
#endif

void wf::Vectors::begin_cycle(uint_t icycle) {
    reduction::clear_local(m_stats.m_summed);
    m_store.attempt_remap();
    m_refs.begin_cycle(icycle);
}

void wf::Vectors::end_cycle(uint_t icycle) {
    reduction::all_sum(m_stats.m_summed);
    m_refs.end_cycle(icycle);
}

wf_comp_t wf::Vectors::reference_projected_energy(uint_t ipart) const {
    wf_comp_t num = 0.0;
    wf_comp_t den = 0.0;
    auto& row = m_store.m_row;
    for (auto irec: m_irec_ref_conns) {
        row.jump(irec);
        if (!row.m_ref_conn.get(ipart)) continue;
        DEBUG_ASSERT_FALSE(row.is_freed(), "reference-connected row should not be freed");
        const wf_t& weight = row.m_weight[ipart];
        num += m_ham.get_element(m_refs[ipart].mbf(), row.m_mbf) * weight;
        if (row.m_mbf == m_refs[ipart].mbf()) den+=weight;
    }
    num = mpi::all_sum(num);
    den = mpi::all_sum(den);
    DEBUG_ASSERT_NE(std::abs(den), 0.0, "reference weight is zero");
    return num / den;
}

wf_comp_t wf::Vectors::debug_square_norm(uint_t ipart) const {
    wf_comp_t res = 0.0;
    auto fn = [&](const Walker& row) {
        const wf_t& weight = row.m_weight[ipart];
        res += std::pow(std::abs(weight), 2.0);
    };
    m_store.foreach_row_in_use(fn);
    return mpi::all_sum(res);
}

wf_comp_t wf::Vectors::debug_l1_norm(uint_t ipart) const {
    wf_comp_t res = 0.0;
    auto fn = [&](const Walker& row) {
        const wf_t& weight = row.m_weight[ipart];
        res += std::abs(weight);
    };
    m_store.foreach_row_in_use(fn);
    return mpi::all_sum(res);
}

void wf::Vectors::set_weight(Walker& walker, uint_t ipart, wf_t new_weight) {
    DEBUG_ASSERT_FALSE(std::isnan(std::abs(new_weight)), "new weight is invalid");
    if (m_ref_weights_preserved && walker.m_mbf==m_refs[ipart].mbf()) return;
    wf_t& weight = walker.m_weight[ipart];
    m_stats.m_nwalker.delta()[ipart] += std::abs(new_weight) - std::abs(weight);
    m_stats.m_l2_norm_square.delta()[ipart] += std::pow(std::abs(new_weight), 2.0) - std::pow(std::abs(weight), 2.0);
    weight = new_weight;
}

void wf::Vectors::change_weight(Walker& walker, uint_t ipart, wf_t delta) {
    set_weight(walker, ipart, walker.m_weight[ipart] + delta);
}

void wf::Vectors::scale_weight(Walker& walker, uint_t ipart, double factor) {
    set_weight(walker, ipart, factor * walker.m_weight[ipart]);
}

void wf::Vectors::zero_weight(Walker& walker, uint_t ipart) {
    set_weight(walker, ipart, 0.0);
}

void wf::Vectors::remove_row(Walker& walker) {
    DEBUG_ASSERT_TRUE(m_store.lookup(walker.m_mbf), "MBF doesn't exist in table!");
    for (uint_t ipart = 0ul; ipart < m_format.m_nelement; ++ipart) {
        zero_weight(walker, ipart);
        --m_stats.m_nocc_mbf.delta();
    }
    remove_ref_conn(walker);
    m_store.erase(walker.m_mbf);
}

void wf::Vectors::add_ref_conn(const Walker& walker) {
    if (walker.m_ref_conn.is_clear()) return;
    DEBUG_ASSERT_FALSE(m_irec_ref_conns.count(walker.index()), "this record index is already in the set");
    m_irec_ref_conns.insert(walker.index());
}

void wf::Vectors::remove_ref_conn(const Walker& walker) {
    if (walker.m_ref_conn.is_clear()) return;
    DEBUG_ASSERT_TRUE(m_irec_ref_conns.count(walker.index()), "this record index should be in the set");
    m_irec_ref_conns.erase(walker.index());
}

Walker& wf::Vectors::create_row_(uint_t icycle, const Mbf& mbf, tag::Int<1>) {

    DEBUG_ASSERT_TRUE(mpi::i_am(m_dist.irank(mbf)),
                      "this method should only be called on the rank responsible for storing the MBF");
    auto& row = m_store.insert(mbf);
    ++m_stats.m_nocc_mbf.delta();
    DEBUG_ASSERT_EQ(row.key_field(), mbf, "MBF was not properly copied into key field of WF row");
    row.m_hdiag = m_ham.get_energy(mbf);
    /*
     * we need to be very careful here of off-by-one-like mistakes. the initial walker is "created" at the beginning
     * of MC cycle 0, and so the stats line output for cycle 0 will show that the number of walkers is the initial
     * occupation of the initial row. if a spawning event leads to the creation of another row, it is created on
     * iteration 1 even though it is added in the annihilating call of iteration 0. so, if this method is called in
     * the annihilating process of MC cycle i, it actually "becomes occupied" on cycle i+1.
     */
    if (storing_av_weights()) {
        row.m_icycle_occ = icycle+1;
        row.m_average_weight = 0;
    }
    return row;
}

Walker& wf::Vectors::create_row_(uint_t icycle, const Mbf& mbf, tag::Int<0>) {
    auto& row = create_row_(icycle, mbf, tag::Int<1>());
    for (uint_t ipart=0ul; ipart < npart(); ++ipart) {
        row.m_ref_conn.put(ipart, m_refs[ipart].connected(mbf));
    }
    add_ref_conn(row);
    return row;
}

Spawn& wf::Vectors::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator,
                              bool deterministic, uint_t dst_ipart) {
    auto& dst_table = send(m_dist.irank(dst_mbf));

    auto& spawn = dst_table.m_row;
    spawn.push_back_jump();

    spawn.m_dst_mbf = dst_mbf;
    spawn.m_delta_weight = delta;
    spawn.m_src_initiator = initiator;
    spawn.m_src_deterministic = deterministic;
    spawn.m_ipart_dst = dst_ipart;
    return spawn;
}

Spawn& wf::Vectors::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator, bool deterministic,
                              uint_t dst_ipart, const field::Mbf& src_mbf, wf_t src_weight) {
    auto& spawn = add_spawn(dst_mbf, delta, initiator, deterministic, dst_ipart);
    if (spawn.m_send_parents) {
        spawn.m_src_mbf = src_mbf;
        spawn.m_src_weight = src_weight;
    }
    DEBUG_ASSERT_NE(dst_mbf, src_mbf, "spawning diagonally");
    return spawn;
}

void wf::Vectors::refresh_all_hdiags() {
    auto fn = [&](Walker& row) {
        row.m_hdiag = m_ham.get_energy(row.m_mbf);
    };
    m_store.foreach_row_in_use(fn);
}

void wf::Vectors::refresh_all_ref_conns() {
    m_irec_ref_conns.clear();
    auto fn = [&](Walker& row) {
        for (uint_t ipart=0ul; ipart < npart(); ++ipart) {
            auto connected = m_refs[ipart].connected(row.m_mbf);
            row.m_ref_conn.put(ipart, connected);
        }
        add_ref_conn(row);
    };
    m_store.foreach_row_in_use(fn);
}

void wf::Vectors::fci_init(FciInitOptions opts, uint_t max_ncomm) {
    /*
     * perform the eigensolver procedure for the required number of states
     */
    FciInitializer init(m_ham, opts);
    const auto results = init.solve();
    /*
     * compute the ratio of initial number of walkers to L1-norms of the eigenvectors to get the right scale
     */
    v_t<wf_t> scale_facs;
    if (mpi::i_am_root()) {
        v_t<ham_t> evals;
        results.get_evals(evals);
        logging::info("FCI energies ({} root{}): {}", nroot(), string::plural(nroot()), convert::to_string(evals));

        const auto nw = m_opts.m_wavefunction.m_nw_init.m_value;
        for (uint_t iroot=0ul; iroot<opts.m_nroot; ++iroot)
            scale_facs.push_back(nw / math::l1_norm(results.get_evec(iroot), results.nelement_evec()));
    }

    uint_t irow = 0ul;
    /*
     * continue to distribute the eigenvectors in blocks until there are no remaining elements
     */
    char done = false;
    while (!mpi::all_land(done)) {
        if (mpi::i_am_root()) {
            auto& row = init.m_mbf_order_table.m_row;
            auto& mbf = row.m_field;
            const auto irow_end = std::min(init.m_mbf_order_table.nrow_in_use(), irow + max_ncomm);
            for (row.jump(irow); row.in_range(irow_end); ++row) {
                for (uint_t iroot = 0ul; iroot < nroot(); ++iroot) {
                    for (uint_t ireplica = 0ul; ireplica < nreplica(); ++ireplica) {
                        auto ipart = m_format.flatten({iroot, ireplica});
                        const auto weight = results.get_evec(iroot)[row.index()]*scale_facs[iroot];
                        add_spawn(mbf, weight, true, false, ipart);
                    }
                }
            }
            irow = irow_end;
            done = (irow == init.m_mbf_order_table.nrow_in_use());
        } else {
            done = true;
        }
        /*
         * use the spawning send/recv tables to distribute the wavefunction from the root rank to the correct ranks
         */
        m_send_recv.communicate();
        auto& recv_row = m_send_recv.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            m_store.lookup(recv_row.m_dst_mbf, m_store.m_insert_row);
            if (!m_store.m_insert_row) create_row(0ul, recv_row.m_dst_mbf);
            m_store.m_insert_row.m_weight = recv_row.m_delta_weight;
        }
    }
}

void wf::Vectors::orthogonalize(reduction::NdArray<wf_t, 3>& overlaps, uint_t iroot, uint_t jroot, uint_t ireplica) {
    ASSERT(iroot <= jroot);
    auto& row = m_store.m_row;
    const auto ipart_src = m_format.flatten({iroot, ireplica});
    const auto ipart_dst = m_format.flatten({jroot, ireplica});
    overlaps.m_local[{iroot, jroot, ireplica}] +=
            arith::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_dst];
    if (jroot + 1 < nroot()) {
        // there is another part to project onto
        const auto ipart_next = m_format.flatten({jroot + 1, ireplica});
        overlaps.m_local[{iroot, jroot + 1, ireplica}] +=
                arith::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_next];
    }
    if (iroot < jroot) {
        const auto& overlap = overlaps.m_reduced[{iroot, jroot, ireplica}];
        const auto& norm = overlaps.m_reduced[{iroot, iroot, ireplica}];
        ASSERT(std::abs(norm) > 1e-12);
        const auto gs_coeff = overlap / norm;
        change_weight(row, ipart_dst, -gs_coeff * row.m_weight[ipart_src]);
    }
}

void wf::Vectors::orthogonalize() {
    // bra root, ket root, replica
    reduction::NdArray<wf_t, 3> overlaps({nroot(), nroot(), nreplica()});
    auto& row = m_store.m_row;
    for (uint_t iroot = 0ul; iroot < nroot(); ++iroot) {
        for (uint_t jroot = iroot; jroot < nroot(); ++jroot) {
            for (uint_t ireplica = 0ul; ireplica < nreplica(); ++ireplica) {
                for (row.restart(); row; ++row) {
                    if (!row.m_mbf.is_clear()) orthogonalize(overlaps, iroot, jroot, ireplica);
                }
                overlaps.all_sum();
            }
        }
    }
}

void wf::Vectors::save(const hdf5::NodeWriter& parent) const {
    auto& row = m_store.m_row;
    hdf5::GroupWriter gw(parent, "wf");
    row.m_mbf.save(gw, true);
    row.m_weight.save(gw, true);
}

void wf::Vectors::save() const {
    REQUIRE_TRUE_ALL(m_opts.m_wavefunction.m_save.m_enabled, "wavefunction saving is disabled in config document")
    hdf5::FileWriter fw(m_opts.m_wavefunction.m_save.m_path);
    save(fw);
}
