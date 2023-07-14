//
// Created by Robert John Anderson on 2020-04-03.
//

#include "M7_lib/basis/Suites.h"
#include "Wavefunction.h"
#include "CiInitializer.h"
#include "HfExcitHists.h"
#include "M7_lib/hdf5/DistTableLoader.h"


v_t<TableBase::Loc> wf::Vectors::setup() {
    v_t<TableBase::Loc> ref_locs;
    /*
     * create the reference MBF and add it to the walker table
     */
    buffered::Mbf ref_mbf(m_sector);
    mbf::set(ref_mbf, m_sector.particles(), m_opts.m_reference.m_mbf_init, 0ul);


    const auto pmntr = m_opts.m_wavefunction.m_ci_pmntr.m_enabled;

    if (pmntr) {
        /*
         * using the "permanitiator" adaptation with low-rank CI information as a source
         */
        hf_excit_hist::initialize(*this, ref_mbf, m_opts.m_wavefunction.m_ci_pmntr);
    }

    /*
     * if the input-specified nw is 0, assume we will initialize the WF with the target number of walkers
     */
    const wf_t nw_init = (m_opts.m_wavefunction.m_nw_init.m_value == 0.0) ?
            m_opts.m_shift.m_nw_targets.m_value[0] : m_opts.m_wavefunction.m_nw_init.m_value;
    /*
     * insert reference MBF into the store table
     */
    {
        const auto ref_loc = create_row_setup(0, ref_mbf);
        if (ref_loc.is_mine()) {
            auto ref_walker = m_store.m_row;
            ref_walker.jump(ref_loc.m_irec);
            for (uint_t ipart = 0ul; ipart < npart(); ++ipart) set_weight(ref_walker, ipart, 1.0);
        }

        for (auto ipart=0ul; ipart<npart(); ++ipart) ref_locs.push_back(ref_loc);
    }

    const auto& init_space_type = m_opts.m_wavefunction.m_init_space.m_type.m_value;
    const auto init_space_solve = m_opts.m_wavefunction.m_init_space.m_solve.m_value;

    if (m_opts.m_wavefunction.m_load_large_ci.m_enabled) {
        // the wavefunction (MBFs only) is to be loaded from HDF5 archive
        hdf5::FileReader fr(m_opts.m_wavefunction.m_load_large_ci.m_path);
        load(fr);
    }
    if (m_opts.m_wavefunction.m_load.m_enabled) {
        // the wavefunction is to be loaded from HDF5 archive
        hdf5::FileReader fr(m_opts.m_wavefunction.m_load.m_path);
        load(fr);
    }
    else if (init_space_type != "ref") {
        // the wavefunction is to be initialized using eigenvectors from the Arnoldi method
        logging::info("Performing exact CI initialization of wavefunctions");
    }
    {
        ci_init::Options opts;
        opts.m_nroot = this->nroot();
        if (init_space_type == "fci") {
            opts.m_loop_kind = ci_init::Options::Conns;
            ci_init::FciSubspace subspace(&m_ham, m_sector.particles());
            ci_init(subspace, opts, init_space_solve);
        }
        else if (init_space_type == "ref_conn") {
            opts.m_loop_kind = ci_init::Options::MbfPairs;
            ci_init::RefConnSubspace subspace(&m_ham, ref_mbf);
            ci_init(subspace, opts, init_space_solve);
        }
        else if (init_space_type == "ref") {
            auto flipped = ref_mbf;
            flipped.ms2_flip();
            const auto flip_fac = m_opts.m_wavefunction.m_init_space.m_ms2_flip;
            if (flip_fac && flipped != ref_mbf) {
                const auto flipped_loc = create_row_setup(0, flipped);
                if (flipped_loc.is_mine()) {
                    auto flipped_walker = m_store.m_row;
                    flipped_walker.jump(flipped_loc.m_irec);
                    for (uint_t ipart = 0ul; ipart < npart(); ++ipart) set_weight(flipped_walker, ipart, flip_fac);
                }
            }
        }
    }

    /*
     * scale each part to the required initial number of walkers
     */
    for (uint_t ipart = 0ul; ipart < npart(); ++ipart) {
        const auto& ref_loc = ref_locs[ipart];
        auto scale_fac = nw_init;
        if (m_opts.m_shift.m_fix_ref_weight) {
            // initial and target walker numbers pertain to the reference population
            wf_t ref_weight = 0.0;
            if (ref_loc.is_mine()) {
                m_store.m_row.jump(ref_loc.m_irec);
                ref_weight = m_store.m_row.m_weight[ipart];
            }
            mpi::bcast(ref_weight, ref_loc.m_irank);
            REQUIRE_TRUE_ALL(ref_weight, "reference should have been found with non-zero weight");
            scale_fac /= ref_weight;
        }
        else {
            // initial and target walker numbers pertain to the total population
            scale_fac /= debug_l1_norm(ipart);
        }
        auto fn = [&](Walker& row) {scale_weight(row, ipart,  scale_fac);};
        m_store.foreach_row_in_use(fn);
    }

    return ref_locs;
}
v_t<double> wf::Vectors::make_stoch_thresh_mags() const {
    const auto& input = m_opts.m_wavefunction.m_stoch_thresh_mags.m_value;
    if (input.empty()) return v_t<double>(nshift_space(), 1.0);
    REQUIRE_EQ_ALL(input.size(), nshift_space(),
                   "if stochastic threshold magnitudes are specified at input, there must be one per shift space");
    auto any_neg = std::any_of(input.cbegin(), input.cend(), [](double v){return v<0;});
    REQUIRE_FALSE_ALL(any_neg, "all specified stochastic threshold magnitudes should be non-negative");
    return input;
}

uint_t wf::Vectors::nshift_space() const {
    return m_opts.m_shift.m_nw_targets.m_value.size();
}

wf::Vectors::Vectors(const conf::Document& opts, const Hamiltonian& ham):
    communicator::BasicSend<Walker, Spawn>(
        "wavefunction",
        // walker row:
        {
            ham.m_basis,
            opts.m_wavefunction.m_nroot,
            opts.m_av_ests.need_replication() ? 2ul:1ul
        },
        opts.m_wavefunction.m_distribution,
        // store sizing
        {
            uint_t(opts.m_shift.nw_target_total()),
            opts.m_wavefunction.m_buffers.m_store_exp_fac
        },
        // send/recv row
        {ham.m_basis, need_send_parents(opts)},
        // send/recv sizing
        {
            std::max(10ul, uint_t(opts.m_shift.nw_target_total() * opts.m_propagator.m_tau_init)),
            opts.m_wavefunction.m_buffers.m_comm_exp_fac
        }
    ),
    m_opts(opts),
    m_ham(ham),
    m_sector(m_ham.m_basis, m_ham.default_particles(m_opts.m_particles)),
    m_format(m_store.m_row.m_weight.m_format),
    m_stats(m_format, nshift_space()),
    m_large_ci_set(m_opts.m_wavefunction.m_large_ci_set.m_enabled ?
        new mbf::table_t("large CI set", mbf::row_t({m_ham.m_basis, Walker::c_mbf_field_name})) : nullptr),
    m_gathered_hist(MbfWeightRow(m_store.m_row), false),
    m_stoch_round_mags(make_stoch_thresh_mags()),
    m_refs(opts.m_reference, *this, setup()),
    m_chkpt_files(opts.m_wavefunction.m_chkpt){

    REQUIRE_TRUE(m_send_recv.recv().m_row.m_dst_mbf.belongs_to_row(), "row-field reference error");

    logging::info("Distributing wavefunction rows in {} block{}", m_dist.nblock(),
                  string::plural(m_dist.nblock()));
    if (m_large_ci_set) {
        logging::info("Keeping list of all MBFs which at any point remain occupied for >= {} with average weight >= {}",
                      string::plural("cycle", m_opts.m_wavefunction.m_large_ci_set.m_ncycle_thresh),
                      m_opts.m_wavefunction.m_large_ci_set.m_av_weight_thresh);
        m_large_ci_set->set_expansion_factor(m_store.get_expansion_factor());
    }
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

    const auto l2_norm_square = m_stats.m_l2_norm_square.total()[ipart];
    REQUIRE_GT(l2_norm_square, 0.0, "L2 norm must be positive non-zero");
    for (row.restart(); row; ++row) {
        rows.push_back({
            std::to_string(row.index()),
            row.m_mbf.to_string(),
            convert::to_string(row.m_weight[ipart], {true, 6}),
            convert::to_string(row.m_weight[ipart] / std::sqrt(l2_norm_square), {false, 4}),
            convert::to_string(row.exceeds_initiator_thresh(ipart, m_opts.m_propagator.m_nadd)),
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
    if (m_large_ci_set) {
        auto& row = m_large_ci_set->m_row;
        hdf5::FileWriter fw(m_opts.m_wavefunction.m_large_ci_set.m_path);
        hdf5::GroupWriter gw(fw, "wf");
        row.m_field.save(gw, true);
    }
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

void wf::Vectors::begin_cycle(uint_t icycle) {
    reduction::clear_local(m_stats.m_summed);
    m_store.remap_if_due();
    m_refs.begin_cycle(icycle);
}

void wf::Vectors::end_cycle(uint_t icycle) {
    attempt_chkpt(icycle);
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

uint_t wf::Vectors::debug_ndeterministic(uint_t iroot) const {
    uint_t res = 0;
    auto fn = [&](const Walker& row) {res += row.m_deterministic.get(iroot);};
    m_store.foreach_row_in_use(fn);
    return res;
}

void wf::Vectors::set_weight(Walker& walker, uint_t ipart, wf_t new_weight, uint_t new_shift_space) {
    DEBUG_ASSERT_FALSE(math::is_nan_or_inf(std::abs(new_weight)), "new weight is invalid");
    if (m_ref_weights_preserved && walker.m_mbf==m_refs[ipart].mbf()) return;
    wf_t& weight = walker.m_weight[ipart];
    const auto delta = std::abs(new_weight) - std::abs(weight);
    m_stats.m_nw.delta()[ipart] += delta;
    {
        // update number of walkers resolved by shift space index
        const uint_t old_shift_space = walker.m_shift_space;
        const auto& format = m_stats.m_nw_by_shift_space.m_format;
        /*
         * ipart is a compound index of root and replica, so combine with the minor index to obtain the overall flat
         * index for the shift space
         */
        auto iflat_old = format.combine<1>(old_shift_space, ipart);
        auto iflat_new = format.combine<1>(new_shift_space, ipart);
        if (iflat_old == iflat_new) m_stats.m_nw_by_shift_space.delta()[iflat_new] += delta;
        else {
            m_stats.m_nw_by_shift_space.delta()[iflat_old] -= std::abs(weight);
            m_stats.m_nw_by_shift_space.delta()[iflat_new] += std::abs(new_weight);
            --m_stats.m_nocc_mbf_by_shift_space.delta()[old_shift_space];
            ++m_stats.m_nocc_mbf_by_shift_space.delta()[new_shift_space];
            walker.m_shift_space = new_shift_space;
            // protect if the new space is S0, there exist higher spaces, and the walker is a reference connection
            const auto nspace = m_stats.m_nocc_mbf_by_shift_space.m_format.m_nelement;
            if (new_shift_space==0 && (nspace > 1) && walker.m_ref_conn.get(ipart)) walker.protect();
        }
    }
    m_stats.m_l2_norm_square.delta()[ipart] += std::pow(std::abs(new_weight), 2.0) - std::pow(std::abs(weight), 2.0);
    weight = new_weight;
}

void wf::Vectors::change_weight(Walker& walker, uint_t ipart, wf_t delta, uint_t new_shift_space) {
    set_weight(walker, ipart, walker.m_weight[ipart] + delta, new_shift_space);
}

void wf::Vectors::scale_weight(Walker& walker, uint_t ipart, double factor, uint_t new_shift_space) {
    set_weight(walker, ipart, factor * walker.m_weight[ipart], new_shift_space);
}

void wf::Vectors::scale_weight(Walker& walker, uint_t ipart, double factor) {
    scale_weight(walker, ipart, factor, walker.m_shift_space);
}

void wf::Vectors::zero_weight(Walker& walker, uint_t ipart) {
    set_weight(walker, ipart, 0.0);
}

void wf::Vectors::remove_row(Walker& walker) {
    DEBUG_ASSERT_TRUE(m_store.lookup(walker.m_mbf), "MBF doesn't exist in table!");
    for (uint_t ipart = 0ul; ipart < m_format.m_nelement; ++ipart) {
        zero_weight(walker, ipart);
        --m_stats.m_nocc_mbf.delta();
        --m_stats.m_nocc_mbf_by_shift_space.delta()[walker.m_shift_space];
    }
    remove_ref_conn(walker);
    m_store.erase(walker.m_mbf);
}

void wf::Vectors::try_add_to_large_ci_set(Walker& walker, uint_t icycle, Epochs& shift_epoch) {
    if (!m_large_ci_set) return;
    if (icycle < shift_epoch.icycle_start_last() + m_opts.m_wavefunction.m_large_ci_set.m_delay) return;
    const auto occ_ncycle = walker.occupied_ncycle(icycle);
    if (occ_ncycle < m_opts.m_wavefunction.m_large_ci_set.m_ncycle_thresh.m_value) return;
    const auto av_weight = walker.m_average_weight[0] / occ_ncycle;
    if (std::abs(av_weight) < m_opts.m_wavefunction.m_large_ci_set.m_av_weight_thresh.m_value) return;
    if (m_large_ci_set->lookup(walker.m_mbf)) return;
    m_large_ci_set->insert(walker.m_mbf);
    ++m_stats.m_nlarge_ci.delta();
}

void wf::Vectors::discretize(Walker& walker, PRNG& prng) {
    const auto round_mag = m_stoch_round_mags[walker.m_shift_space];
    // no rounding to be done if the thresh is 0
    if (round_mag == 0.0) return;
    // leave weight alone if the walker is protected from deletion
    if (walker.is_protected()) return;
    for (uint_t ipart=0ul; ipart < walker.m_wf_format.m_nelement; ++ipart) {
        // retrieve the post-death weight
        const auto weight = walker.m_weight[ipart];
        // don't attempt stochastic round if the weight exceeds the threshold
        if (std::abs(weight) >= round_mag) return;
        // else, do the stochastic round, logging the change in magnitude
        const auto new_weight = prng.stochastic_round(weight, round_mag);
        set_weight(walker, ipart, new_weight);
    }
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

Walker& wf::Vectors::create_row_(uint_t icycle, const Mbf& mbf, uint_t shift_space, tag::Int<1>) {
    DEBUG_ASSERT_TRUE(mpi::i_am(m_dist.irank(mbf)),
                      "this method should only be called on the rank responsible for storing the MBF");
    auto& row = m_store.insert(mbf);
    ++m_stats.m_nocc_mbf.delta();
    ++m_stats.m_nocc_mbf_by_shift_space.delta()[shift_space];
    DEBUG_ASSERT_EQ(row.key_field(), mbf, "MBF was not properly copied into key field of WF row");
    row.m_hdiag = m_ham.get_energy(mbf);
    row.m_shift_space = shift_space;
    // protect if the new space is S0 and there exist higher spaces
    const auto nspace = m_stats.m_nocc_mbf_by_shift_space.m_format.m_nelement;
    if (shift_space==0 && (nspace > 1)) row.protect();
    row.m_log_enhancement_fac = 0.0;
    /*
     * we need to be very careful here of off-by-one-like mistakes. the initial walker is "created" at the beginning
     * of MC cycle 0, and so the stats line output for cycle 0 will show that the number of walkers is the initial
     * occupation of the initial row. if a spawning event leads to the creation of another row, it is created on
     * iteration 1 even though it is added in the annihilating call of iteration 0. so, if this method is called in
     * the annihilating process of MC cycle i, it actually "becomes occupied" on cycle i+1.
     */
    row.m_icycle_occ = icycle+1;
    row.m_average_weight = 0;
    return row;
}

Walker& wf::Vectors::create_row_(uint_t icycle, const Mbf& mbf, uint_t shift_space, tag::Int<0>) {
    if (m_opts.m_wavefunction.m_no_row_creation) {
        m_store.m_row.select_null();
        return m_store.m_row;
    }
    auto& row = create_row_(icycle, mbf, shift_space, tag::Int<1>());
    for (uint_t ipart=0ul; ipart < npart(); ++ipart) {
        auto is_ref_conn = m_refs[ipart].connected(mbf);
        row.m_ref_conn.put(ipart, is_ref_conn);
        // all reference connections are automatically in shift space 0
        const auto nspace = m_stats.m_nocc_mbf_by_shift_space.m_format.m_nelement;
        if (is_ref_conn && (nspace > 1)) change_weight(row, ipart, 0.0, 0);
    }
    add_ref_conn(row);
    return row;
}

Spawn& wf::Vectors::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator,
                              bool deterministic, uint_t dst_ipart, uint_t dst_shift_space) {
    auto& dst_table = send(m_dist.irank(dst_mbf));

    auto& spawn = dst_table.m_row;
    spawn.push_back_jump();

    spawn.m_dst_mbf = dst_mbf;
    spawn.m_delta_weight = delta;
    spawn.m_src_initiator = initiator;
    spawn.m_src_deterministic = deterministic;
    spawn.m_ipart_dst = dst_ipart;
    spawn.m_dst_shift_space = dst_shift_space;
    return spawn;
}

Spawn& wf::Vectors::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator, bool deterministic,
                              uint_t dst_ipart, const field::Mbf& src_mbf, wf_t src_weight, uint_t dst_shift_space) {
    auto& spawn = add_spawn(dst_mbf, delta, initiator, deterministic, dst_ipart, dst_shift_space);
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

void wf::Vectors::ci_init(const ci_init::Subspace& subspace, v_t<const wf_t*> weight_vecs, int ms2_flip_fac, uint_t max_ncomm) {
    char have_weights = !weight_vecs.empty();
    mpi::bcast(have_weights);

    const auto& table = subspace.m_mbf_order_table;
    buffered::Mbf flipped(m_sector);

    uint_t irow = 0ul;
    /*
     * continue to distribute the eigenvectors in blocks until there are no remaining elements
     */
    char done = false;
    while (!mpi::all_land(done)) {
        if (mpi::i_am_root()) {
            auto& row = table.m_row;
            auto& mbf = row.m_field;
            const auto irow_end = std::min(table.nrow_in_use(), irow + max_ncomm);
            for (row.jump(irow); row.in_range(irow_end); ++row) {
                for (uint_t iroot = 0ul; iroot < nroot(); ++iroot) {
                    const auto weight_vec = weight_vecs.empty() ? nullptr : weight_vecs[iroot];
                    for (uint_t ireplica = 0ul; ireplica < nreplica(); ++ireplica) {
                        auto ipart = m_format.flatten({iroot, ireplica});
                        const auto weight = weight_vec ? weight_vec[row.index()] : 0.0;
                        if (!have_weights || std::abs(weight) > 1e-6) {
                            add_spawn(mbf, weight, true, false, ipart, 0);
                            if (ms2_flip_fac) {
                                flipped = mbf;
                                mbf::ms2_flip(flipped);
                                if (flipped != mbf) add_spawn(flipped, ms2_flip_fac * weight, true, false, ipart, 0);
                            }
                        }
                    }
                }
            }
            irow = irow_end;
            done = (irow == table.nrow_in_use());
        } else {
            done = true;
        }
        /*
         * use the spawning send/recv tables to distribute the wavefunction from the root rank to the correct ranks
         */
        m_send_recv.communicate();
        auto& recv_row = m_send_recv.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            auto& store_row = lookup_or_create_row_setup_(0, recv_row.m_dst_mbf);
            if (have_weights) store_row.m_weight = recv_row.m_delta_weight;
        }
    }
}

void wf::Vectors::ci_init(const ci_init::Subspace& subspace, ci_init::Options opts, uint_t max_ncomm) {
    /*
     * perform the eigensolver procedure for the required number of states
     */
    ci_init::Initializer init(subspace, opts);
    const auto results = init.solve();

    if (mpi::i_am_root()) {
        v_t<ham_t> evals;
        results.get_evals(evals);
        logging::info("CI energies ({} root{}): {}", nroot(), string::plural(nroot()), convert::to_string(evals));
    }
    ci_init(subspace, results.get_evecs(), m_opts.m_wavefunction.m_init_space.m_ms2_flip, max_ncomm);
}

void wf::Vectors::ci_init(const ci_init::Subspace& subspace, uint_t max_ncomm) {
    ci_init(subspace, v_t<const wf_t*>(), m_opts.m_wavefunction.m_init_space.m_ms2_flip, max_ncomm);
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

void wf::Vectors::load(const hdf5::NodeReader& parent) {
    hdf5::GroupReader gr(parent, "wf");
    uintv_t weight_shape = this->m_store.m_row.m_weight.m_format.shape_vector();
    const uint_t nreplica = this->nreplica();
    auto nreplica_on_file = nreplica;
    const auto have_weights = gr.child_exists("weight");

    if (have_weights) {
        weight_shape = hdf5::DatasetLoader::read_format(gr, "weight", true, true).m_local.m_item.m_shape;
        nreplica_on_file = weight_shape.back();

        REQUIRE_EQ_ALL(weight_shape.front(), nroot(), "incompatible number of roots in file");

        if (nreplica > nreplica_on_file) {
            logging::info("Loading non-replicated wavefunctions for a replica calculation: duplicating weights");
        } else if (nreplica < nreplica_on_file) {
            logging::warn("Loading replicated wavefunctions for a non-replica calculation: discarding second replica");
        }
    }
    else logging::info("Loading file with only MBFs specified, no weights");

    auto file_ipart_fn = [&nreplica, &nreplica_on_file](uint_t ipart) {
        if (nreplica == nreplica_on_file) return ipart;
        else if (nreplica > nreplica_on_file) return ipart / 2;
        return ipart * 2;
    };

    struct LoadRow : Row {
        NdFormat<c_ndim_wf> m_format;
        field::Mbf m_mbf;
        field::Numbers<wf_t, c_ndim_wf> m_weight;
        LoadRow(sys::Basis basis, uintv_t weight_shape):
            m_format(array::from_vector<uint_t, c_ndim_wf>(weight_shape)),
            m_mbf(this, basis, Walker::c_mbf_field_name),
            m_weight(this, m_format, "weight"){}
    };

    typedef buffered::Table<LoadRow> load_table_t;
    load_table_t load_table("WF load table", {m_sector.basis(), weight_shape});

    v_t<FieldBase*> fields = {&load_table.m_row.m_mbf};
    if (have_weights) fields.push_back(&load_table.m_row.m_weight);
    DistTableLoader loader(gr, fields);

    // total number of rows received in all communications
    uint_t nrow_recv = 0ul;

    auto fill_fn = [&](uint_t nitem) {
        auto& row = load_table.m_row;
        for (row.restart(); row.in_range(nitem); ++row) {
            auto irank_dst = m_dist.irank(row.m_mbf);
            auto& send_table = send(irank_dst);
            auto& send_row = send_table.m_row;
            for (uint_t ipart=0ul; ipart < npart(); ++ipart) {
                send_row.push_back_jump();
                send_row.m_dst_mbf = row.m_mbf;
                send_row.m_ipart_dst = ipart;
                const auto file_ipart = file_ipart_fn(ipart);
                if (have_weights) send_row.m_delta_weight = row.m_weight[file_ipart];
            }
        }
        m_send_recv.communicate();

        auto fn = [&](const Spawn& recv_row) {
            auto& store_row = lookup_or_create_row_setup_(0, recv_row.m_dst_mbf);
            store_row.protect();
            const auto ipart = recv_row.m_ipart_dst[0];
            if (have_weights) set_weight(store_row, ipart, recv_row.m_delta_weight);
            ++nrow_recv;
        };
        recv().foreach_row_in_use(fn);
    };
    const uint_t nitem_per_op = 100000;
    logging::info("Loading walkers from HDF5 archive (upto {} items per read operation)", nitem_per_op);
    logging::info_("Reading {} items locally, {} items globally", loader.nitem_local(), loader.nitem());
    loader.load(nitem_per_op, fill_fn);
    REQUIRE_EQ_ALL(mpi::all_sum(nrow_recv), loader.nitem(), "not all walkers loaded");
    logging::info("{} wavefunction rows successfully loaded from HDF5 archive", loader.nitem());
    logging::info("{} total wavefunction rows", mpi::all_sum(m_store.nrow_in_use()));
}

bool wf::Vectors::was_loaded() const {
    return m_opts.m_wavefunction.m_load.m_enabled;
}

void wf::Vectors::attempt_chkpt(uint_t icycle) {
    if (!m_chkpt_files) return;
    auto path = m_chkpt_files.get_file_path(icycle);
    if (path.empty()) return;
    logging::info("Saving wavefunction to checkpoint file {} on cycle {}", path, icycle);
    hdf5::FileWriter fw(path);
    fw.save_attr("icycle", icycle);
    save(fw);
}

void wf::Vectors::update_gathered_hist(wf_comp_t thresh, uint_t icycle) {
    logging::info("Gathering histogrammed CI weights");
    if (thresh == 0.0) logging::info("Not discarding based on average weight");
    else logging::info("Discarding MBFs with average weight < {} from histogrammed set", thresh);
    logging::flush_all();

    uint_t ndiscard = 0ul;
    buffered::Table<MbfWeightRow> local_averaged(MbfWeightRow{m_store.m_row});
    auto& local_row = local_averaged.m_row;
    local_row.restart();

    auto add_local_hist_walker_fn = [&ndiscard, &thresh, &local_row, &icycle](const Walker& walker){
        if (walker.is_protected()) {
            const auto av_weight = walker.m_average_weight[0] / walker.occupied_ncycle(icycle);
            if (std::abs(av_weight) < thresh) {
                ++ndiscard;
                return;
            }
            local_row.push_back_jump();
            local_row.m_mbf = walker.m_mbf;
            local_row.m_weight = walker.m_average_weight;
        }
    };
    m_store.foreach_row_in_use(add_local_hist_walker_fn);

    logging::info("Local histogrammed rows collected - performing all MPI all gatherv");
    logging::flush_all();
    // TODO: node-shared gathered_averaged
    m_gathered_hist.all_gatherv(local_averaged);

    ndiscard = mpi::all_sum(ndiscard);
    if (ndiscard) logging::info("Discarded {} low-weight MBFs from the histogrammed set", ndiscard);
    logging::flush_all();

    m_last_gathered_hist_thresh = thresh;
    m_last_gathered_hist_icycle = icycle;
}

void wf::Vectors::update_gathered_hist_if_changed(wf_comp_t thresh, uint_t icycle) {
    if (m_gathered_hist.empty() ||
        thresh != m_last_gathered_hist_thresh ||
        icycle != m_last_gathered_hist_icycle)
        update_gathered_hist(thresh, icycle);
}

void wf::Vectors::attempt_gathered_hist_save(uint_t icycle) {
    if (!m_opts.m_wavefunction.m_save_hist.m_enabled) return;
    update_gathered_hist_if_changed(m_opts.m_wavefunction.m_save_hist.m_thresh, icycle);
    hdf5::FileWriter fw(m_opts.m_wavefunction.m_save_hist.m_path);
    auto& row = m_gathered_hist.m_row;
    hdf5::GroupWriter gw(fw, "wf");
    row.m_mbf.save(gw, mpi::i_am_root());
    row.m_weight.save(gw, mpi::i_am_root());
}
