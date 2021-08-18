//
// Created by rja on 10/11/2020.
//

#include "Solver.h"

Solver::Solver(const fciqmc_config::Document &opts, Propagator &prop, Wavefunction &wf,
               std::vector<TableBase::Loc> ref_locs) :
        m_prop(prop), m_opts(prop.m_opts), m_wf(wf),
        m_refs(m_opts.m_reference, m_prop.m_ham, m_wf, ref_locs),
        m_exit("exit"),
        m_uniform_twf(
                m_opts.m_inst_ests.m_spf_uniform_twf ?
                new UniformTwf(m_prop.m_ham, m_wf.npart(), prop.m_ham.nsite()) : nullptr),
        m_weighted_twf(m_opts.m_inst_ests.m_spf_weighted_twf ?
                       new WeightedTwf(m_prop.m_ham, m_wf.npart(), prop.m_ham.nsite(),
                                       m_opts.m_inst_ests.m_spf_weighted_twf.m_fermion_fac,
                                       m_opts.m_inst_ests.m_spf_weighted_twf.m_boson_fac) : nullptr),
        m_maes(m_opts.m_av_ests, m_prop.m_ham.nsite(), m_prop.m_ham.nelec()),
        m_annihilator(m_wf, m_prop, m_refs, m_maes.m_bilinears.m_rdms, m_icycle, opts.m_propagator.m_nadd),
        m_archive(opts) {

    log::info("Replicating walker populations: {}", m_wf.nreplica()==2);
    if (m_wf.nreplica()==2 && !m_prop.nexcit_gen())
        log::warn("Replica populations are redundant when doing exact propagation");

    if (m_maes.m_bilinears && m_wf.nreplica() == 1 && m_prop.nexcit_gen())
        log::warn("Attempting a stochastic propagation estimation of bilinear MAEs without replication, "
                  "this is biased");

    if (mpi::i_am_root()) {
        m_stats = std::unique_ptr<FciqmcStats>(new FciqmcStats("M7.stats", "FCIQMC", {m_prop}));
        m_timing_stats = std::unique_ptr<TimingStats>(new TimingStats("M7.timing", "FCIQMC Timings", {}));
    }
    if (m_opts.m_stats.m_parallel)
        m_parallel_stats = std::unique_ptr<ParallelStats>(
                new ParallelStats("M7.stats." + std::to_string(mpi::irank()), "FCIQMC Parallelization", {}));

    /**
     * setup archive members
     */
    m_archive.add_member(m_prop);
    if (m_maes.m_ref_excits) m_archive.add_member(m_maes.m_ref_excits);
    if (m_maes.m_bilinears) {
        if (m_maes.m_bilinears.m_rdms) m_archive.add_member(m_maes.m_bilinears.m_rdms);
        if (m_maes.m_bilinears.m_spec_moms) m_archive.add_member(m_maes.m_bilinears.m_spec_moms);
    }
    /**
     * read previous calculation data into archivable objects if archive loading is enabled
     */
    m_archive.load();

    m_wf.m_ra.activate(m_icycle);
}

void Solver::execute(size_t ncycle) {
    log::info("Beginning solver loop...");
    for (size_t i = 0ul; i < ncycle; ++i) {
        m_cycle_timer.reset();
        m_cycle_timer.unpause();
        begin_cycle();

        m_propagate_timer.reset();
        m_propagate_timer.unpause();
        if (m_detsub) m_detsub->update();
        loop_over_occupied_mbfs();
        m_propagate_timer.pause();

        m_communicate_timer.reset();
        m_communicate_timer.unpause();

        m_wf.communicate();
        m_communicate_timer.pause();

        m_annihilate_timer.reset();
        m_annihilate_timer.unpause();
        loop_over_spawned();
        mpi::barrier();
        if (m_detsub) {
            m_detsub->make_rdm_contribs(m_maes.m_bilinears.m_rdms, m_refs[0].get_mbf());
            m_detsub->project(m_prop.tau());
        }
        m_annihilate_timer.pause();

        end_cycle();
        m_cycle_timer.pause();
        output_stats();

        m_maes.output(m_icycle, m_prop.m_ham);
        ++m_icycle;

        if (m_exit.read() && m_exit.m_v) {
            log::info("exit requested from file, terminating solver loop at MC cycle {}", i);
            break;
        }
        if (m_maes.m_accum_epoch) {
            if (i == m_maes.m_accum_epoch.icycle_start() + m_opts.m_av_ests.m_ncycle) {
                if (m_icycle == ncycle)
                    log::info("maximum number of MEV accumulating cycles ({}) "
                              "reached at MC cycle {}", m_opts.m_av_ests.m_ncycle, i);
                break;
            }
        }
        log::flush();
    }
    if (m_icycle == ncycle) log::info("maximum cycle number ({}) reached", m_icycle);
    if (m_maes.m_accum_epoch) {
        // repeat the last cycle but do not perform any propagation
        --m_icycle;
        finalizing_loop_over_occupied_mbfs();
    }
    std::cout << m_maes.m_bilinears.m_rdms.get_energy(m_prop.m_ham.m_frm)+7.878453652277122 << std::endl;
}

void Solver::begin_cycle() {
    m_chk_nwalker_local = m_wf.m_nwalker.m_local[{0, 0}] + m_wf.m_delta_nwalker.m_local[{0, 0}];
    m_chk_ninitiator_local = m_wf.m_ninitiator.m_local[{0, 0}] + m_wf.m_delta_ninitiator.m_local[{(0, 0)}];
    m_wf.begin_cycle();
    ASSERT(m_wf.m_nwalker.m_local[0] == 0);
    m_wf.m_ra.update(m_icycle);
    m_propagate_timer.reset();
    if (m_wf.nroot() > 1 && m_prop.m_shift.m_variable_mode) {
        m_wf.orthogonalize();
    }
    m_refs.begin_cycle();

    auto update_epoch = [&](const size_t &ncycle_wait) {
        if (m_wf.nreplica() == 2) {
            if (m_prop.m_shift.m_variable_mode[0] && m_prop.m_shift.m_variable_mode[1]) {
                auto max_start = std::max(
                        m_prop.m_shift.m_variable_mode[0].icycle_start(),
                        m_prop.m_shift.m_variable_mode[1].icycle_start());
                if (m_icycle > max_start + ncycle_wait) return true;
            }
        } else {
            if (m_prop.m_shift.m_variable_mode[0]) {
                auto start = m_prop.m_shift.m_variable_mode[0].icycle_start();
                if (m_icycle > start + ncycle_wait) return true;
            }
        }
        return false;
    };

    if (m_maes) {
        if (m_maes.m_accum_epoch.update(m_icycle, update_epoch(m_opts.m_av_ests.m_delay))) {
            REQUIRE_TRUE_ALL(m_maes.all_stores_empty(),
                             "MAEs only beginning to be accumulated, but not all store tables are empty");
        }
    }

    if (m_opts.m_propagator.m_semistochastic && !m_detsub) {
        if (update_epoch(m_opts.m_propagator.m_semistochastic.m_delay)) {
            m_detsub = std::unique_ptr<DeterministicSubspace>(new DeterministicSubspace(
                    m_opts.m_propagator.m_semistochastic, m_wf, m_icycle));
            m_detsub->build_from_most_occupied(m_prop.m_ham, m_maes.m_bilinears);
            //m_detsub->build_from_all_occupied(m_prop.m_ham);
            log::debug("Initialized deterministic subspace");
        }
    }
}

void Solver::loop_over_occupied_mbfs() {
    auto &row = m_wf.m_store.m_row;

    for (row.restart(); row.in_range(); row.step()) {
        /*
         * stats always refer to the state of the wavefunction in the previous iteration
         */

        if (row.m_mbf.is_zero())
            /*
             * this is a free row, caused by an earlier call to m_wf.remove_row()
             */
            continue;


        if (row.m_weight.is_zero() && !row.is_protected()) {
            /*
             * MBF has become unoccupied in all parts and must be removed from mapped list, but it must first make all
             * associated averaged contributions to MEVs
             */
            m_maes.make_average_contribs(row, m_refs, m_icycle);
            m_wf.remove_row();
            continue;
        }

        /*
         * if the accumulation of MEVs has just started, treat the row as though it became occupied in the annihilation
         * loop of the last MC cycle.
         */
        if (m_maes.m_accum_epoch.started_this_cycle(m_icycle)) {
            DEBUG_ASSERT_TRUE(m_maes.m_accum_epoch, "should be in MAE accumulation epoch");
            row.m_icycle_occ = m_icycle - 1;
            row.m_average_weight = 0;
        }

        if (m_maes.m_accum_epoch) {
            row.m_average_weight += row.m_weight;
        }

        if (m_maes.is_period_cycle(m_icycle)) {
            /*
             * this is the end of a planned block-averaging cycle, therefore there may be unaccounted-for contributions
             * which need to be included in the average
             */
            m_maes.make_average_contribs(row, m_refs, m_icycle);
        }

        m_refs.contrib_row();

        for (size_t ipart = 0ul; ipart < m_wf.m_format.m_nelement; ++ipart) {

            DEBUG_ASSERT_TRUE(!m_wf.m_store.m_row.m_mbf.is_zero(),
                              "Stored MBF should not be zeroed");
            DEBUG_ASSERT_TRUE(mpi::i_am(m_wf.get_rank(m_wf.m_store.m_row.m_mbf)),
                              "Stored MBF should be on its allocated rank");

            const auto &weight = row.m_weight[ipart];

            m_wf.m_nocc_mbf.m_local++;
            if (row.m_initiator.get(ipart))
                m_wf.m_ninitiator.m_local[ipart]++;

            m_wf.m_nwalker.m_local[ipart] += std::abs(weight);
            m_wf.m_l2_norm_square.m_local[ipart] += std::pow(std::abs(weight), 2.0);

            if (ipart == 0) {
                if (m_opts.m_inst_ests.m_spf_uniform_twf) m_uniform_twf->add(row.m_weight, row.m_mbf);
                if (m_opts.m_inst_ests.m_spf_weighted_twf) m_weighted_twf->add(row.m_weight, row.m_mbf);
            }

            if (m_wf.m_ra.is_active()) {
                m_spawning_timer.reset();
                m_spawning_timer.unpause();
            }
            propagate_row(ipart);
            if (m_wf.m_ra.is_active()) {
                m_spawning_timer.pause();
                m_wf.m_ra.record_work_time(row, m_spawning_timer);
            }
        }
    }
    m_synchronization_timer.reset();
    m_synchronization_timer.unpause();
    mpi::barrier();
    m_synchronization_timer.pause();
}

void Solver::finalizing_loop_over_occupied_mbfs() {
    if (!m_maes.m_accum_epoch) return;
    auto &row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        if (!row.m_mbf.is_zero()) m_maes.make_average_contribs(row, m_refs, m_icycle);
    }
    m_maes.end_cycle();
}

void Solver::loop_over_spawned() {
    mpi::barrier();
    if (!m_wf.recv().m_hwm) {
        log::debug_("no rows received, omitting loop over spawned");
        return;
    }
    auto hwm_before = m_wf.m_store.m_hwm;
    DEBUG_ONLY(hwm_before);

    m_annihilator.sort_recv();
    m_annihilator.loop_over_dst_mbfs();

    DEBUG_ASSERT_LE(m_wf.m_store.m_hwm, hwm_before + m_wf.recv().m_hwm,
                    "the store table shouldn't have grown by more rows than were received!");
    m_wf.recv().clear();
}

void Solver::propagate_row(const size_t &ipart) {
    auto &row = m_wf.m_store.m_row;

    if (row.is_cleared()) return;

    if (consts::float_is_zero(row.m_weight[ipart])) return;

    m_prop.off_diagonal(m_wf, ipart);
    m_prop.diagonal(m_wf, ipart);
}

void Solver::end_cycle() {
    /*
     * TODO: make these checks compatible with dynamic rank allocation
     */
//    double chk_ratio;
    if (!consts::float_is_zero(m_wf.m_nwalker.m_local[{0, 0}])) {
//        chk_ratio = m_chk_nwalker_local / m_wf.m_nwalker(0, 0);
//        bool chk = m_chk_nwalker_local == 0.0 || consts::floats_nearly_equal(chk_ratio, 1.0);
//        if (!chk) std::cout << "discrepancy: " << m_chk_nwalker_local-m_wf.m_nwalker(0, 0) << std::endl;
//        MPI_REQUIRE(chk,"Unlogged walker population changes have occurred");
    }
//    MPI_REQUIRE(m_chk_ninitiator_local == m_wf.m_ninitiator(0, 0),
//                "Unlogged creations of initiator MBFs have occurred");
    m_refs.end_cycle();
    //if (m_mevs.m_fermion_rdm) m_mevs.m_fermion_rdm->end_cycle();
    m_wf.end_cycle();
    m_prop.update(m_icycle, m_wf);
    if (m_uniform_twf) m_uniform_twf->reduce();
    if (m_weighted_twf) m_weighted_twf->reduce();
}

void Solver::output_stats() {

    auto sync_overhead = mpi::all_sum(static_cast<double>(m_synchronization_timer));
    if (mpi::i_am_root()) {
        auto &stats = m_stats->m_row;
        stats.m_icycle = m_icycle;
        stats.m_tau = m_prop.tau();
        stats.m_shift = m_prop.m_shift.m_values;
        stats.m_nwalker = m_wf.m_nwalker.m_reduced;
        stats.m_delta_nwalker = m_wf.m_delta_nwalker.m_reduced;
        stats.m_nwalker_spawned = m_wf.m_nspawned.m_reduced;
        stats.m_nwalker_annihilated = m_wf.m_nannihilated.m_reduced;
        stats.m_ref_proj_energy_num = m_refs.proj_energy_nums();
        stats.m_ref_weight = m_refs.weights();
        stats.m_ref_proj_energy = stats.m_ref_proj_energy_num;
        stats.m_ref_proj_energy /= stats.m_ref_weight;
        stats.m_l2_norm = m_wf.m_l2_norm_square.m_reduced;
        stats.m_l2_norm.to_sqrt();
        stats.m_ninitiator = m_wf.m_ninitiator.m_reduced;
        stats.m_nocc_mbf = m_wf.m_nocc_mbf.m_reduced;
        stats.m_delta_nocc_mbf = m_wf.m_delta_nocc_mbf.m_reduced;
        if (m_prop.nexcit_gen()) stats.m_exlvl_probs = m_prop.exlvl_probs();
        if (m_uniform_twf) stats.m_uniform_twf_num = m_uniform_twf->m_numerator_total[0];
        if (m_weighted_twf) {
            stats.m_weighted_twf_num = m_weighted_twf->m_numerator_total[0];
            stats.m_weighted_twf_denom = m_weighted_twf->m_denominator_total[0];
        }
        stats.m_reweighting_factor = m_prop.m_shift.m_reweighter.m_total;
        m_stats->flush();

        auto &timing_stats = m_timing_stats->m_row;
        timing_stats.m_total_synchronization_overhead = sync_overhead;
        timing_stats.m_propagate_loop_time = m_propagate_timer;
        timing_stats.m_communication_time = m_communicate_timer;
        timing_stats.m_annihilation_loop_time = m_annihilate_timer;
        timing_stats.m_total_cycle_time = m_cycle_timer;
        m_timing_stats->flush();
    }

    if (m_opts.m_stats.m_parallel) {
        auto &stats = m_parallel_stats->m_row;
        stats.m_icycle = m_icycle;
        stats.m_synchronization_overhead = m_synchronization_timer;
        stats.m_nblock_wf_ra = m_wf.m_ra.nblock_local();
        stats.m_nwalker_total = m_wf.m_nwalker.m_reduced.sum();
        stats.m_nwalker_lookup_skip = m_wf.m_store.m_nskip_total;
        stats.m_nwalker_lookup = m_wf.m_store.m_nlookup_total;
        stats.m_nrow_recv = m_wf.m_comm.m_last_recv_count;
        m_parallel_stats->flush();
    }
}