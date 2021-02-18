//
// Created by rja on 10/11/2020.
//

#include "Solver.h"

void Solver::loop_over_occupied_onvs() {
    /*
     * Loop over all rows in the m_wf.m_walkers table and if the row is not empty:
     *      ascertain the initiator status of the ONV
     *      if not initiator and weight > initiator threshold, then grant initiator status
     *      if initiator and weight < initiator threshold, then revoke initiator status
     *      perform all the off-diagonal propagation (fill send table)
     *      update local weight in the diagonal cloning/death step
     * else if all elements of the m_weight field are zero, the row should be removed
     */
    auto& row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()){
        /*
         * stats always refer to the state of the wavefunction in the previous iteration
         */

        if (row.is_cleared()) {
            // row is empty
            continue;
        }

        MPI_ASSERT(!m_wf.m_store.m_row.m_onv.is_zero(),
                   "Stored ONV should not be zeroed");
        MPI_ASSERT(mpi::i_am(m_wf.get_rank(m_wf.m_store.m_row.m_onv)),
                   "Stored ONV should be on its allocated rank");

        const auto weight = row.m_weight(m_wf.m_ipart);

        if (consts::float_nearly_zero(weight, 1e-12)){
            // ONV has become unoccupied and must be removed from mapped list
            m_wf.remove_walker();
            continue;
        }

        m_wf.m_nocc_onv(0, 0)++;
        if (row.m_initiator.get(m_wf.m_ipart))
            m_wf.m_ninitiator(0, 0)++;
        m_wf.m_nwalker(0, 0) += std::abs(weight);
        m_wf.m_l2_norm_square(0, 0) += std::pow(std::abs(weight), 2.0);

        m_reference.add_row();

        if (m_wf.m_ra.is_active()) {
            m_spawning_timer.reset();
            m_spawning_timer.unpause();
        }
        propagate_row();
        if (m_wf.m_ra.is_active()) {
            m_spawning_timer.pause();
            m_wf.m_ra.record_work_time(row, m_spawning_timer);
        }
    }
    m_synchronization_timer.reset();
    m_synchronization_timer.unpause();
    mpi::barrier();
    m_synchronization_timer.pause();
}

void Solver::annihilate_row() {
    const auto& row = m_wf.recv().m_row;
    const auto& dst_onv = row.m_dst_onv;
    ASSERT(!dst_onv.is_zero());
    // check that the received determinant has come to the right place
    ASSERT(m_wf.m_ra.get_rank(dst_onv) == mpi::irank())
    const defs::wf_t& delta_weight = row.m_delta_weight;
    // zero magnitude weights should not have been communicated
    ASSERT(!consts::float_is_zero(delta_weight));

    auto irow_walkers = *m_wf.m_store[dst_onv];

#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "bitstring:             " << dst_onv.to_string() << std::endl;
    std::cout << consts::verb << "delta weight:          " << delta_weight << std::endl;
    std::cout << consts::verb << "found in walker list:  " << string_utils::yn(irow_walkers!=~0ul) << std::endl;
#endif

    if (irow_walkers == ~0ul) {
        /*
         * the destination determinant is not currently occupied, so initiator rules
         * must be applied
         */
        if (!row.m_src_initiator) {
            //m_aborted_weight += std::abs(*delta_weight);
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs <<
                "ABORTED SPAWN TO UNOCCUPIED DETERMINANT: PARENT NON-INITIATOR" << std::endl;
#endif
            return;
        }

        m_wf.create_walker_(
                dst_onv,
                delta_weight,
                m_prop.m_ham.get_energy(dst_onv),
                m_reference.is_connected(dst_onv));
    } else {
        m_wf.m_store.m_row.jump(irow_walkers);
        m_wf.change_weight(delta_weight);
    }

//#ifdef VERBOSE_DEBUGGING
//    std::cout << consts::verb << "row in walker list:    " << irow_main << std::endl;
//#endif
//    /*
//     * if we have stochastically generated a connection between determinants in a deterministic
//     * subspace, so we must reject this connection.
//     */
//    if (recv.m_flags.m_src_deterministic(irow_recv) && m_wf.m_walkers.m_flags.m_deterministic(irow_walkers)) {
//#ifdef VERBOSE_DEBUGGING
//        std::cout << consts::verb << consts::chevs <<
//                  "ABORTED SPAWN: STOCHASTICALLY GENERATED A DETERMINISTIC CONNECTION" << std::endl;
//#endif
//        return;
//    }
//    auto weight = m_wf.m_walkers.m_weight(irow_walkers, 0, 0);
//    //m_square_norm.m_delta += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
//    //m_nwalker.m_delta -= std::abs(*weight);
//    weight += delta_weight;
//    //m_nwalker.m_delta += std::abs(*weight);
}

void Solver::loop_over_spawned() {
    mpi::barrier();
    const auto &row = m_wf.recv().m_row;
    for (row.restart(); row.in_range(); row.step()){
        annihilate_row();
    }
    m_wf.recv().clear();
}

Solver::Solver(Propagator &prop, Wavefunction &wf, Table::Loc ref_loc) :
        m_prop(prop),
        m_opts(prop.m_opts),
        m_wf(wf),
        m_reference(m_opts, m_prop.m_ham, m_wf, 0, ref_loc) {
    if (mpi::i_am_root())
        m_stats = mem_utils::make_unique<StatsFile<FciqmcStatsSpecifier>>("M7.stats");
    m_parallel_stats = mem_utils::make_unique<StatsFile<ParallelStatsSpecifier>>(
            "M7.stats."+std::to_string(mpi::irank()));
}

void Solver::execute(size_t niter) {
    for (size_t i = 0ul; i < niter; ++i) {
        m_cycle_timer.reset();
        m_cycle_timer.unpause();
        begin_cycle();

        m_propagate_timer.reset();
        m_propagate_timer.unpause();
        loop_over_occupied_onvs();
        m_propagate_timer.pause();

        m_communicate_timer.reset();
        m_communicate_timer.unpause();
        m_wf.communicate();
        m_communicate_timer.pause();

        m_annihilate_timer.reset();
        m_annihilate_timer.unpause();
        loop_over_spawned();
        m_annihilate_timer.pause();

        end_cycle();
        m_cycle_timer.pause();
        output_stats();
        ++m_icycle;
    }
}

void Solver::propagate_row() {
    auto& row = m_wf.m_store.m_row;
    const auto& ipart = m_wf.m_ipart;

    if (row.is_cleared()) return;

    //const auto onv = m_wf.m_walkers.m_onv(irow);
    const auto weight = row.m_weight(ipart);
    if (consts::float_is_zero(weight)) return;

//    bool is_initiator = m_wf.m_walkers.m_flags.m_initiator(irow, 0, 0);
    bool is_deterministic = row.m_deterministic.get(ipart);
//    bool is_ref_connection = m_wf.m_walkers.m_flags.m_reference_connection(irow);

//    if (consts::float_is_zero(weight) && !is_deterministic) {
//#ifdef VERBOSE_DEBUGGING
//        std::cout << consts::verb << consts::chevs << "ZERO WEIGHT: REMOVING FROM LIST" << std::endl;
//                std::cout << consts::verb << "is initiator:     " << flag_initiator << std::endl;
//                std::cout << consts::verb << "weight:           " << *weight << std::endl;
//#endif
//        if (is_initiator) {
//#ifdef VERBOSE_DEBUGGING
//            std::cout << consts::verb << consts::chevs << "INITIATOR STATUS REVOKED: DETERMINANT REMOVED" << std::endl;
//#endif
//            //m_ninitiator.m_delta--;
//        }
//        //m_nocc_det.m_delta--;
//        //m_data.remove(irow);
//        return;
//    }
//
//
//    if (!is_initiator && std::abs(weight) >= m_opts.nadd_initiator) {
//#ifdef VERBOSE_DEBUGGING
//        std::cout << consts::verb << consts::chevs << "INITIATOR STATUS GRANTED" << std::endl;
//#endif
//        is_initiator = true;
//        //m_ninitiator.m_delta++;
//    }
//    /*
//    else if (flag_initiator && std::abs(*weight) < m_input.nadd_initiator) {
//#ifdef VERBOSE_DEBUGGING
//            std::cout << consts::verb << consts::chevs << "INITIATOR STATUS REVOKED: WEIGHT FELL BELOW THRESHOLD MAGNITUDE" << std::endl;
//#endif
//        // initiator status revoked
//        // flag_initiator = false;
//        // delta_ninitiator--;
//    }
//     */


    m_prop.off_diagonal(m_wf);
    m_prop.diagonal(m_wf);

    if (!is_deterministic && consts::float_is_zero(weight)) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "ALL WALKERS DIED: REMOVING DETERMINANT FROM LIST" << std::endl;
#endif
//                if (flag_initiator) m_ninitiator.m_delta--;
//                m_nocc_det.m_delta--;
//                m_data.remove(irow);
    }
}

void Solver::begin_cycle() {
    m_chk_nwalker_local = m_wf.m_nwalker(0, 0) + m_wf.m_delta_nwalker(0, 0);
    m_chk_ninitiator_local = m_wf.m_ninitiator(0, 0) + m_wf.m_delta_ninitiator(0, 0);
    m_wf.begin_cycle();
    if (m_prop.m_variable_shift.started_last_cycle(m_icycle))
        m_wf.m_ra.activate(m_icycle);
    m_wf.m_ra.update(m_icycle);
    m_propagate_timer.reset();
    m_reference.begin_cycle();
}

void Solver::end_cycle() {
    //double chk_ratio;
    if (!consts::float_is_zero(m_wf.m_nwalker(0, 0))) {
        //chk_ratio = m_chk_nwalker_local / m_wf.m_nwalker(0, 0);
        //if (m_chk_nwalker_local > 0.0 && !consts::floats_nearly_equal(chk_ratio, 1.0))
        //    throw std::runtime_error("Unlogged walker population changes have occurred");
    }

    if (m_chk_ninitiator_local != m_wf.m_ninitiator(0, 0)) {
        //throw std::runtime_error("Unlogged creations of initiator ONVs have occurred");
    }

    m_wf.end_cycle();
    m_reference.end_cycle();
    m_prop.update(m_icycle, m_wf);
}

void Solver::output_stats() {

    auto sync_overhead = mpi::all_sum((double)m_synchronization_timer);
    if (mpi::i_am_root()) {
        m_stats->m_icycle() = m_icycle;
        m_stats->m_tau() = m_prop.tau();
        m_stats->m_shift() = m_prop.m_shift;
        m_stats->m_nwalker() = m_wf.m_nwalker.reduced(0, 0);
        m_stats->m_delta_nwalker() = m_wf.m_delta_nwalker.reduced(0, 0);
        m_stats->m_ref_proj_energy_num() = m_reference.proj_energy_num();
        m_stats->m_ref_weight() = m_reference.m_weight(0);
        m_stats->m_ref_proj_energy() = m_reference.proj_energy();
        m_stats->m_l2_norm() = std::sqrt(m_wf.m_l2_norm_square.reduced(0, 0));
        m_stats->m_ninitiator() = m_wf.m_ninitiator.reduced(0, 0);
        m_stats->m_nocc_onv() = m_wf.m_nocc_onv.reduced(0, 0);
        m_stats->m_psingle() = m_prop.m_magnitude_logger.m_psingle;
        m_stats->m_total_synchronization_overhead() = sync_overhead;
        m_stats->m_propagate_loop_time() = m_propagate_timer;
        m_stats->m_communication_time() = m_communicate_timer;
        m_stats->m_annihilation_loop_time() = m_annihilate_timer;
        m_stats->m_total_cycle_time() = m_cycle_timer;
        m_stats->flush();
    }


    m_parallel_stats->m_icycle() = m_icycle;
    m_parallel_stats->m_synchronization_overhead() = m_synchronization_timer;
    m_parallel_stats->m_nblock_wf_ra() = m_wf.m_ra.nblock_();
    m_parallel_stats->m_nwalker() = m_wf.m_nwalker(0, 0);
//    m_parallel_stats->m_nrow_free_walker_list() = m_wf.m_walkers.
//    StatsColumn<size_t> m_walker_list_high_water_mark;
//    StatsColumn<double> m_walker_list_high_water_mark_fraction;
//    StatsColumn<size_t> m_nrow_sent;
//    StatsColumn<size_t> m_largest_nrow_sent;
//    StatsColumn<double> m_largest_send_list_filled_fraction;
//    StatsColumn<size_t> m_irank_largest_nrow_sent;
//    StatsColumn<size_t> m_nrow_recv;
//    StatsColumn<double> m_recv_list_filled_fraction;
    m_parallel_stats->flush();
}
