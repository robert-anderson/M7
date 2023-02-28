//
// Created by Robert John Anderson on 2020-04-03.
//

#include <M7_lib/basis/Suites.h>
#include "Wavefunction.h"
#include "M7_lib/linalg/FciIters.h"
#include "M7_lib/sort/QuickSort.h"
#include "FciInitializer.h"
#include "M7_lib/util/Math.h"

/*
MappedTableBase::nbucket_guess(
        opts.m_propagator.m_nw_target / mpi::nrank(),
        opts.m_wavefunction.m_hash_mapping.m_remap_ratio),
opts.m_wavefunction.m_hash_mapping.m_remap_nlookup,
opts.m_wavefunction.m_hash_mapping.m_remap_ratio
 */

/*
 *
    Communicator(str_t name,
                 const store_row_t &store_row, DistribOptions dist_opts, Sizing store_sizing,
                 const send_recv_row_t &send_recv_row, Sizing comm_sizing):
 */


wf::Fci::Fci(const conf::Document& opts, const sys::Sector& sector) :
    communicator::BasicSend<Walker, Spawn>(
        "wavefunction",
        // walker row:
        {
            sector.basis(),
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
        {sector.basis(), need_send_parents(opts)},
        // send/recv sizing
        {
            std::max(10ul, uint_t(opts.m_propagator.m_nw_target * opts.m_propagator.m_tau_init)),
            opts.m_wavefunction.m_buffers.m_comm_exp_fac
        }
    ),
    m_opts(opts),
    m_sector(sector),
    m_format(m_store.m_row.m_weight.m_format),
    m_ninitiator(m_format),
    m_nwalker(m_format),
    m_delta_nwalker(m_format),
    m_l2_norm_square(m_format),
    m_delta_l2_norm_square(m_format),
    m_nspawned(m_format),
    m_nannihilated(m_format) {

    REQUIRE_TRUE(m_send_recv.recv().m_row.m_dst_mbf.belongs_to_row(), "row-field reference error");
    m_summables.add_members(m_ninitiator, m_nocc_mbf, m_delta_nocc_mbf,
                            m_nwalker, m_delta_nwalker, m_l2_norm_square, m_delta_l2_norm_square,
                            m_nspawned, m_nannihilated);

    logging::info("Distributing wavefunction rows in {} block{}", m_dist.nblock(),
                  string::plural(m_dist.nblock()));
}

void wf::Fci::log_top_weighted(uint_t ipart, uint_t nrow) {
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
    rows.push_back({"", "many-body basis function", "walker number", "initiator", "semistoch", "MPI rank"});
    for (row.restart(); row; ++row) {
        rows.push_back({
            std::to_string(row.index()),
            row.m_mbf.to_string(),
            convert::to_string(row.m_weight[ipart], 6),
            convert::to_string(row.is_initiator(ipart, m_opts.m_propagator.m_nadd)),
            convert::to_string(bool(row.m_deterministic[iroot_part(ipart)])),
            convert::to_string(m_dist.irank(row.m_mbf))
        });
    }
    logging::info_table("Top-weighted WF elements for part "+std::to_string(ipart), rows, true, false, 1ul);
}

wf::Fci::~Fci() {
    for (uint_t ipart=0ul; ipart<npart(); ++ipart) log_top_weighted(ipart);
}

void wf::Fci::h5_write(const hdf5::NodeWriter& parent, str_t name) {
    auto field_names = []() -> strv_t {
        if (c_enable_fermions != c_enable_bosons) return {"mbf", "weight"};
        else if (c_enable_fermions != c_enable_bosons) return {"mbf (fermion)", "mbf (boson)", "weight"};
        return {};
    };
    m_store.save(parent, name, field_names(), true);
}

void wf::Fci::h5_read(const hdf5::NodeReader& /*parent*/, const Hamiltonian& /*ham*/, const field::Mbf& /*ref*/,
                      str_t /*name*/) {
#if 0
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
#endif
}

void wf::Fci::begin_cycle() {
    m_summables.zero_all_local();
    m_store.attempt_remap();
}

void wf::Fci::end_cycle() {
    m_summables.all_sum();
}

wf_comp_t wf::Fci::square_norm(uint_t ipart) const {
    wf_comp_t res = 0.0;
    auto& row = m_store.m_row;
    for (row.restart(); row; ++row) {
        const wf_t& weight = row.m_weight[ipart];
        res += std::pow(std::abs(weight), 2.0);
    }
    return mpi::all_sum(res);
}

wf_comp_t wf::Fci::l1_norm(uint_t ipart) const {
    wf_comp_t res = 0.0;
    auto& row = m_store.m_row;
    for (row.restart(); row; ++row) {
        const wf_t& weight = row.m_weight[ipart];
        res += std::abs(weight);
    }
    return mpi::all_sum(res);
}

void wf::Fci::set_weight(Walker& walker, uint_t ipart, wf_t new_weight) {
    wf_t& weight = walker.m_weight[ipart];
    m_delta_nwalker.m_local[ipart] += std::abs(new_weight);
    m_delta_nwalker.m_local[ipart] -= std::abs(weight);
    m_delta_l2_norm_square.m_local[ipart] += std::pow(std::abs(new_weight), 2.0);
    m_delta_l2_norm_square.m_local[ipart] -= std::pow(std::abs(weight), 2.0);
    weight = new_weight;
}

void wf::Fci::change_weight(Walker& walker, uint_t ipart, wf_t delta) {
    set_weight(walker, ipart, walker.m_weight[ipart] + delta);
}

void wf::Fci::scale_weight(Walker& walker, uint_t ipart, double factor) {
    set_weight(walker, ipart, factor * walker.m_weight[ipart]);
}

void wf::Fci::zero_weight(Walker& walker, uint_t ipart) {
    set_weight(walker, ipart, 0.0);
}

void wf::Fci::remove_row(Walker& walker) {
    DEBUG_ASSERT_TRUE(m_store.lookup(walker.m_mbf), "MBF doesn't exist in table!");
    for (uint_t ipart = 0ul; ipart < m_format.m_nelement; ++ipart) {
        zero_weight(walker, ipart);
        m_delta_nocc_mbf.m_local--;
    }
    m_store.erase(walker.m_mbf);
}

Walker& wf::Fci::create_row_(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, const v_t<bool>& refconns) {
    DEBUG_ASSERT_EQ(refconns.size(), npart(), "should have as many reference rows as WF parts");
    DEBUG_ASSERT_TRUE(mpi::i_am(m_dist.irank(mbf)),
                      "this method should only be called on the rank responsible for storing the MBF");
    auto& row = m_store.insert(mbf);
    m_delta_nocc_mbf.m_local++;
    DEBUG_ASSERT_EQ(row.key_field(), mbf, "MBF was not properly copied into key field of WF row");
    row.m_hdiag = hdiag;
    for (uint_t ipart=0ul; ipart < npart(); ++ipart)
        row.m_ref_conn.put(ipart, refconns[ipart]);
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

TableBase::Loc wf::Fci::create_row(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, const v_t<bool>& refconns) {
    const uint_t irank = m_dist.irank(mbf);
    uint_t irec;
    if (mpi::i_am(irank)) {
        irec = create_row_(icycle, mbf, hdiag, refconns).index();
    }
    mpi::bcast(irec, irank);
    return {irank, irec};
}

Spawn& wf::Fci::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator,
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

Spawn& wf::Fci::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator, bool deterministic,
                          uint_t dst_ipart, const field::Mbf& src_mbf, wf_t src_weight) {
    auto& spawn = add_spawn(dst_mbf, delta, initiator, deterministic, dst_ipart);
    if (spawn.m_send_parents) {
        spawn.m_src_mbf = src_mbf;
        spawn.m_src_weight = src_weight;
    }
    DEBUG_ASSERT_NE(dst_mbf, src_mbf, "spawning diagonally");
    return spawn;
}

void wf::Fci::fci_init(const Hamiltonian& h, FciInitOptions opts, uint_t max_ncomm) {
    /*
     * perform the eigensolver procedure for the required number of states
     */
    FciInitializer init(h, opts);
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
            if (!m_store.m_insert_row) create_row(0ul, recv_row.m_dst_mbf, h.get_energy(recv_row.m_dst_mbf), 0);
            m_store.m_insert_row.m_weight = recv_row.m_delta_weight;
        }
    }
}