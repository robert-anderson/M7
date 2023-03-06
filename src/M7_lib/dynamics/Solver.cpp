//
// Created by Robert J. Anderson on 10/11/2020.
//

#include "Solver.h"

std::unique_ptr<HartreeFock> Solver::make_hf() const {
    /*
     * even if there is no true HF determinant, we need to treat the initial ref as one for averaged excits or when
     * this behaviour is explicitly requested in the configuration document
     */
    const bool force_hf = m_opts.m_reference.m_assume_hf_like.m_value || m_opts.m_av_ests.m_hf_excits.m_enabled;
    /*
     * assuming that if there exists a Hartree-Fock determinant, it will be the initial reference of the zero-th root
     */
    if (force_hf || m_prop.m_ham.has_brillouin_theorem(m_refs[0].mbf())) {
        const TableBase::Loc loc(m_refs[0].irec());
        return ptr::smart::make_unique<HartreeFock>(m_wf.m_store, loc, m_opts.m_hf_excits);
    }
    return nullptr;
}

Solver::Solver(const conf::Document &opts, Propagator &prop, wf::Fci &wf,
               v_t<TableBase::Loc> ref_locs) :
        m_opts(opts), m_prop(prop), m_wf(wf),
        m_refs(m_opts.m_reference, m_prop.m_ham, m_wf, ref_locs),
        m_hf(make_hf()), m_exit("exit"), m_maes(opts.m_av_ests, m_wf.m_sector, m_wf.nroot()),
        m_inst_ests(m_wf.m_sector, &m_refs, opts.m_inst_ests),
        m_annihilator(m_wf, m_prop, m_refs, m_hf.get(), m_maes.m_rdms, m_icycle, opts.m_propagator.m_nadd),
        m_detsubs(opts.m_propagator.m_semistochastic) {

    logging::info("Replicating walker populations: {}", m_wf.nreplica() == 2);
    if (m_wf.nreplica() == 2 && !m_prop.ncase_excit_gen())
        logging::warn("Replica populations are redundant when doing exact propagation");

    if (m_maes.m_rdms && m_wf.nreplica() == 1 && m_prop.ncase_excit_gen())
        logging::warn("Attempting a stochastic propagation estimation of bilinear MAEs without replication, "
                  "this is biased");

    if (mpi::i_am_root()) {
        m_stats = ptr::smart::make_unique<FciqmcStats>(
            "M7.stats", "FCIQMC", FciqmcStatsRow(m_prop, m_inst_ests), m_opts.m_stats.m_period);
        m_timing_stats = ptr::smart::make_unique<TimingStats>(
            "M7.timing", "FCIQMC Timings", TimingStatsRow(), m_opts.m_stats.m_period);
    }
    if (m_opts.m_stats.m_parallel) {
        m_parallel_stats = ptr::smart::make_unique<ParallelStats>("M7.stats." + std::to_string(mpi::irank()),
                 "FCIQMC Parallelization", ParallelStatsRow(), m_opts.m_stats.m_period);
    }

    if (m_prop.ncase_excit_gen()) {
        logging::info("Global PRNG state checksum of stochastic propagator: {}",
                   integer::to_hex_string(m_prop.checksum()));
        logging::info_("Local PRNG state checksum of stochastic propagator: {}",
                   integer::to_hex_string(m_prop.checksum_()));
    }

    if (m_hf) {
        const auto s = m_hf->mbf().to_string();
        if (m_prop.m_ham.has_brillouin_theorem(m_hf->mbf()))
            logging::info("Brillouin theorem condition satisfied: assuming {} as HF state", s);
        else
            logging::info("Brillouin theorem condition not satisfied, but assuming {} as HF-like state", s);
    }
    else {
        logging::info("Brillouin theorem condition not satisfied: no Hartree-Fock state assumed");
    }

    /**
     * setup archive members
     */
    /*
    m_archive.add_member(m_prop);
    if (m_maes.m_hf_excits) m_archive.add_member(m_maes.m_hf_excits);
    if (m_maes.m_bilinears) {
        if (m_maes.m_bilinears.m_rdms) m_archive.add_member(m_maes.m_bilinears.m_rdms);
        if (m_maes.m_bilinears.m_spec_moms) m_archive.add_member(m_maes.m_bilinears.m_spec_moms);
    }
    */


    if (m_maes.m_rdms) {
        if (m_maes.m_rdms.is_energy_sufficient(m_prop.m_ham))
            logging::info("Specified RDM rank signatures are sufficient for variational energy estimation");
        else
            logging::warn("Specified RDM rank signatures are insufficient for variational energy estimation");
    }

    /**
     * set up the wavefunctions to their FCI values if requested
     */
    if (m_opts.m_wavefunction.m_fci_init) {
        logging::info("Performing exact FCI initialization of wavefunctions");
        FciInitOptions fci_init_opts;
        fci_init_opts.m_nroot = m_wf.nroot();
        m_wf.fci_init(m_prop.m_ham, fci_init_opts);
    }

    if (m_opts.m_wavefunction.m_ci_permanitiator.m_enabled) {
        REQUIRE_TRUE(m_hf.get(), "CI expansion-based permanitiators requires a HF determinant");
        hf_excit_hist::initialize(m_wf, m_hf.get()->mbf(), m_opts.m_wavefunction.m_ci_permanitiator);
        m_wf.refresh_all_hdiags(m_prop.m_ham);
    }

    // TODO: activate load balancing
    //m_wf.m_dist.activate(m_icycle);
}

void Solver::execute(uint_t ncycle) {
    logging::info("Beginning solver loop...");
    for (uint_t i = 0ul; i < ncycle; ++i) {
        m_cycle_timer.reset();
        m_cycle_timer.unpause();
        begin_cycle();

        m_propagate_timer.reset();
        m_propagate_timer.unpause();
        if (m_detsubs) {
            m_detsubs.update();
            m_detsubs.make_rdm_contribs(m_maes.m_rdms, m_hf.get());
            m_detsubs.make_spec_mom_contribs(m_maes.m_spec_moms);
        }
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
        if (m_detsubs) m_detsubs.project(m_prop.tau());
        m_annihilate_timer.pause();

        end_cycle();
        m_cycle_timer.pause();
        output_stats();

        m_maes.output(m_icycle, m_prop.m_ham);
        ++m_icycle;

        if (m_exit.read() && m_exit.m_v) {
            logging::info("exit requested from file, terminating solver loop at MC cycle {}", i);
            break;
        }
        if (m_maes.m_accum_epoch) {
            if (i == m_maes.m_accum_epoch.icycle_start() + m_opts.m_av_ests.m_ncycle) {
                logging::info("Maximum number of MAE accumulating cycles ({}) reached at MC cycle {}. "
                              "Terminating solver loop.",
                              m_opts.m_av_ests.m_ncycle, i);
                break;
            }
        }
        logging::flush();
    }
    if (m_icycle == ncycle) logging::info("maximum cycle number ({}) reached", m_icycle);
    if (m_maes.m_accum_epoch) {
        // repeat the last cycle but do not perform any propagation
        finalizing_loop_over_occupied_mbfs(m_icycle - 1);
    }
}

void Solver::begin_cycle() {
//    m_chk_nwalker_local = m_wf.m_nwalker.m_local[{0, 0}] + m_wf.m_delta_nwalker.m_local[{0, 0}];

    m_prop.update(m_icycle, m_wf, m_refs);
    if (m_prop.m_shift.m_variable_mode && m_opts.m_shift.m_fix_ref_weight) {
        m_wf.m_ref = &m_refs[0].mbf();
        m_wf.m_preserve_ref = true;
    }

    m_wf.begin_cycle();
    DEBUG_ASSERT_TRUE(m_wf.m_nwalker.delta().m_local.is_zero(), "Cyclic reducibles should be zeroed before new cycle");

    // TODO: update load balancing
    //    m_wf.m_ra.update(m_icycle);
    m_propagate_timer.reset();
    if (m_wf.nroot() > 1 && m_prop.m_shift.m_variable_mode) {
        m_wf.orthogonalize();
    }
    m_refs.begin_cycle(m_icycle);
    if (m_hf) m_hf->update();

    auto update_epoch = [&](uint_t ncycle_wait) {
        const auto &epochs = m_prop.m_shift.m_variable_mode;
        if (!epochs) return false;
        return m_icycle > epochs.icycle_start_last() + ncycle_wait;
    };

    if (m_maes) {
        if (m_maes.m_accum_epoch.update(m_icycle, update_epoch(m_opts.m_av_ests.m_delay))) {
            REQUIRE_TRUE_ALL(m_maes.all_stores_empty(),
                             "MAEs only beginning to be accumulated, but not all store tables are empty");
        }
    }
    if (m_hf && m_hf->m_excit_accums) {
        m_hf->m_excit_accums.m_accum_epoch.update(m_icycle, update_epoch(m_opts.m_hf_excits.m_delay));
    }

    m_inst_ests.begin_cycle(m_icycle);

    if (m_opts.m_propagator.m_semistochastic.m_enabled && !m_detsubs) {
        if (update_epoch(m_opts.m_propagator.m_semistochastic.m_delay)) {
            m_detsubs.init(m_prop.m_ham, m_maes, m_wf, m_icycle);
            logging::info("Initialized deterministic subspace");
        }
    }
}

void Solver::loop_over_occupied_mbfs() {
    auto& walker = m_wf.m_store.m_row;

    for (walker.restart(); walker; ++walker) {
        /*
         * stats always refer to the state of the wavefunction in the previous iteration
         */

        if (walker.m_mbf.is_zero())
            /*
             * this is a free row, caused by an earlier call to m_wf.remove_row()
             */
            continue;


        if (walker.m_weight.is_zero() && !walker.is_protected()) {
            /*
             * MBF has become unoccupied in all parts and must be removed from mapped list, but it must first make all
             * associated averaged contributions to MEVs
             */
            m_maes.make_average_contribs(walker, m_hf.get(), m_icycle);
            m_wf.remove_row(walker);
            continue;
        }

        /*
         * if the accumulation of MAEs has just started, treat the row as though it became occupied in the annihilation
         * loop of the last MC cycle.
         */
        if (m_maes.m_accum_epoch.started_this_cycle(m_icycle) && m_wf.storing_av_weights()) {
            DEBUG_ASSERT_TRUE(m_maes.m_accum_epoch, "should be in MAE accumulation epoch");
            walker.m_icycle_occ = m_icycle;
            walker.m_average_weight = 0;
        }

        if (m_maes.m_accum_epoch) {
            if (m_wf.storing_av_weights()) walker.m_average_weight += walker.m_weight;
        }

        if (m_prop.m_shift.m_variable_mode) {
            // todo: give hf excit hists accumulation its own epoch
            if (m_hf && m_hf->m_excit_accums) m_hf->m_excit_accums.add(walker.m_mbf, walker.m_weight[0]);
        }

        if (m_maes.is_period_cycle(m_icycle)) {
            /*
             * this is the end of a planned block-averaging cycle, therefore there may be unaccounted-for contributions
             * which need to be included in the average
             */
            m_maes.make_average_contribs(walker, m_hf.get(), m_icycle);
        }

        m_refs.contrib_row(walker);
        m_inst_ests.make_numerator_contribs(walker);

        ++m_wf.m_nocc_mbf.delta();
        for (uint_t ipart = 0ul; ipart < m_wf.m_format.m_nelement; ++ipart) {

            DEBUG_ASSERT_TRUE(!m_wf.m_store.m_row.m_mbf.is_zero(),
                              "Stored MBF should not be zeroed");
            DEBUG_ASSERT_TRUE(mpi::i_am(m_wf.m_dist.irank(m_wf.m_store.m_row.m_mbf)),
                              "Stored MBF should be on its allocated rank");

            bool initiator = walker.exceeds_initiator_thresh(ipart, m_prop.m_nadd_initiator);
            if (!initiator) {
                initiator = walker.m_permanitiator.get(ipart);
                if (initiator) m_wf.m_ninitiator_perma.m_local[ipart]++;
            }
            if (initiator) m_wf.m_ninitiator.m_local[ipart]++;

            // TODO: load balancing work logging
//            if (m_wf.m_ra.is_active()) {
//                m_spawning_timer.reset();
//                m_spawning_timer.unpause();
//            }
            propagate_row(walker, ipart, initiator);
//            if (m_wf.m_ra.is_active()) {
//                m_spawning_timer.pause();
//                m_wf.m_ra.record_work_time(row, m_spawning_timer);
//            }
        }
    }
    m_synchronization_timer.reset();
    m_synchronization_timer.unpause();
    mpi::barrier();
    m_synchronization_timer.pause();
}

void Solver::finalizing_loop_over_occupied_mbfs(uint_t icycle) {
    if (!m_maes.m_accum_epoch || m_maes.is_period_cycle(icycle)) return;
    auto& walker = m_wf.m_store.m_row;
    for (walker.restart(); walker; ++walker) {
        if (walker.m_mbf.is_zero()) continue;
        m_maes.make_average_contribs(walker, m_hf.get(), icycle);
    }
    m_maes.end_cycle();
    m_maes.output(m_icycle, m_prop.m_ham, true);
}

void Solver::loop_over_spawned() {
    mpi::barrier();
    if (m_wf.recv().empty()) return;

    auto nrow_in_use_before = m_wf.m_store.nrow_in_use();
    DEBUG_ONLY(nrow_in_use_before);

    m_annihilator.sort_recv();
    m_annihilator.loop_over_dst_mbfs();

    DEBUG_ASSERT_LE(m_wf.m_store.nrow_in_use(), nrow_in_use_before + m_wf.recv().nrow_in_use(),
                    "the store table shouldn't have grown by more rows than were received!");
    m_wf.recv().clear();
}

void Solver::propagate_row(Walker& walker, uint_t ipart, bool initiator) {
    if (walker.is_freed()) return;
    if (walker.m_weight[ipart] == 0.0) return;
    m_prop.off_diagonal(m_wf, walker, ipart, initiator);
    m_prop.diagonal(m_wf, walker, ipart);
}

void Solver::end_cycle() {
    /*
     * TODO: make these checks compatible with dynamic rank allocation
     */
//    double chk_ratio;
//    if (!fptol::numeric_zero(m_wf.m_nwalker.m_local[{0, 0}])) {
//        chk_ratio = m_chk_nwalker_local / m_wf.m_nwalker(0, 0);
//        bool chk = m_chk_nwalker_local == 0.0 || dtype::near_equal(chk_ratio, 1.0);
//        if (!chk) std::cout << "discrepancy: " << m_chk_nwalker_local-m_wf.m_nwalker(0, 0) << std::endl;
//        MPI_REQUIRE(chk,"Unlogged walker population changes have occurred");
//    }
//    MPI_REQUIRE(m_chk_ninitiator_local == m_wf.m_ninitiator(0, 0),
//                "Unlogged creations of initiator MBFs have occurred");
    m_refs.end_cycle(m_icycle);
    m_wf.end_cycle();
    m_maes.end_cycle();
    m_inst_ests.end_cycle(m_icycle);
    if (m_hf && m_hf->m_excit_accums) m_hf->m_excit_accums.attempt_chkpt(m_icycle);
}

void Solver::output_stats() {

    auto sync_overhead = mpi::all_sum(static_cast<double>(m_synchronization_timer));
    if (mpi::i_am_root()) {
        auto &stats = m_stats->m_row;
        stats.m_icycle = m_icycle;
        stats.m_tau = m_prop.tau();
        stats.m_shift = m_prop.m_shift.m_values;
        stats.m_nwalker = m_wf.m_nwalker.prev_total();
        stats.m_delta_nwalker = m_wf.m_nwalker.prev_delta();
        stats.m_nwalker_spawned = m_wf.m_nspawned.m_reduced;
        stats.m_nwalker_annihilated = m_wf.m_nannihilated.m_reduced;
        stats.m_ref_proj_energy_num = m_refs.proj_energy_nums();
        stats.m_ref_weight = m_refs.weights();
        stats.m_ref_proj_energy = stats.m_ref_proj_energy_num;
        stats.m_ref_proj_energy /= stats.m_ref_weight;
        stats.m_l2_norm = m_wf.m_l2_norm_square.prev_total();
        stats.m_l2_norm.to_sqrt();
        stats.m_ninitiator = m_wf.m_ninitiator.m_reduced;
        stats.m_ninitiator_perma = m_wf.m_ninitiator_perma.m_reduced;
        stats.m_nocc_mbf = m_wf.m_nocc_mbf.prev_total();
        stats.m_delta_nocc_mbf = m_wf.m_nocc_mbf.prev_delta();
        if (m_prop.ncase_excit_gen()) stats.m_exlvl_probs = m_prop.excit_gen_case_probs();
        if (m_inst_ests.m_spin_square) stats.m_spin_square_num = m_inst_ests.m_spin_square->m_est.m_proj_num.m_reduced;
        m_stats->commit();

        auto &timing_stats = m_timing_stats->m_row;
        timing_stats.m_total_synchronization_overhead = sync_overhead;
        timing_stats.m_propagate_loop_time = m_propagate_timer;
        timing_stats.m_communication_time = m_communicate_timer;
        timing_stats.m_annihilation_loop_time = m_annihilate_timer;
        timing_stats.m_total_cycle_time = m_cycle_timer;
        m_timing_stats->commit();
    }

    if (m_opts.m_stats.m_parallel) {
        auto &stats = m_parallel_stats->m_row;
        stats.m_icycle = m_icycle;
        stats.m_synchronization_overhead = m_synchronization_timer;
        stats.m_nblock_wf_ra = m_wf.m_dist.nblock_();
        stats.m_nwalker_total = m_wf.m_nwalker.prev_total().sum();
        stats.m_nwalker_lookup_skip = m_wf.m_store.m_nskip_total;
        stats.m_nwalker_lookup = m_wf.m_store.m_nlookup_total;
        stats.m_nrow_recv = m_wf.m_send_recv.m_last_recv_count;
        m_parallel_stats->commit();
    }
}
