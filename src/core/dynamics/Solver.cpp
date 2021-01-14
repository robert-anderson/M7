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
    for (size_t irow = 0ul; irow < m_wf.m_walkers.m_hwm; ++irow) {
        /*
         * stats always refer to the state of the wavefunction in the previous iteration
         */

        if (m_wf.m_walkers.is_cleared(irow)) {
            // row is empty
            continue;
        }

        const auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);

        if (consts::float_nearly_zero(weight, 1e-12)){
            // ONV has become unoccupied and must be removed from mapped list
            m_wf.remove_walker(irow);
            continue;
        }

        m_wf.m_nocc_onv(0, 0)++;
        if (m_wf.m_walkers.m_flags.m_initiator(irow, 0, 0))
            m_wf.m_ninitiator(0, 0)++;
        m_wf.m_nwalker(0, 0) += std::abs(weight);
        m_wf.m_l2_norm_square(0, 0) += std::pow(std::abs(weight), 2.0);

        m_reference.add_row(irow);

        propagate_row(irow);
    }
    m_synchronization_timer.reset();
    m_synchronization_timer.unpause();
    mpi::barrier();
    m_synchronization_timer.pause();
}

void Solver::annihilate_row(const size_t &irow_recv) {
    auto &recv = m_wf.m_spawn.recv();
    auto dst_onv = recv.m_dst_onv(irow_recv);
    ASSERT(!dst_onv.is_zero());
    // check that the received determinant has come to the right place
    ASSERT(m_wf.m_ra.get_rank(dst_onv) == mpi::irank())
    auto delta_weight = recv.m_delta_weight(irow_recv, 0, 0);
    // zero magnitude weights should not have been communicated
    ASSERT(!consts::float_is_zero(delta_weight));

    auto irow_walkers = *m_wf.m_walkers[dst_onv];

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
        if (!recv.m_flags.m_src_initiator(irow_recv, 0, 0)) {
            //m_aborted_weight += std::abs(*delta_weight);
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs <<
                "ABORTED SPAWN TO UNOCCUPIED DETERMINANT: PARENT NON-INITIATOR" << std::endl;
#endif
            return;
        }

        m_wf.create_walker(
                dst_onv,
                delta_weight,
                m_prop.m_ham.get_energy(dst_onv),
                m_reference.is_connected(dst_onv));
    } else {
        m_wf.change_weight(irow_walkers, delta_weight);
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
    //m_annihilation_timer.unpause();
    //if (m_prop->semi_stochastic()) m_detsub->update_weights(m_prop->tau(), m_nwalker.m_delta);
    //m_aborted_weight = 0;
    for (size_t irow_recv = 0ul; irow_recv < m_wf.m_spawn.recv().m_hwm; ++irow_recv) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "RECEIVED DETERMINANT" << std::endl;
            std::cout << consts::verb << "row in recv buffer:    " << irow_recv << std::endl;
#endif
        annihilate_row(irow_recv);
    }
    m_wf.m_spawn.recv().clear();
    //mpi::barrier();
    //m_annihilation_timer.pause();
#if 0
    m_data.clear_tombstones();
#endif
}

const Reference &Solver::reference() const {
    return m_reference;
}

Solver::Solver(Propagator &prop, Wavefunction &wf, views::Onv<> ref_onv) :
        m_prop(prop),
        m_opts(prop.m_opts),
        m_wf(wf),
        m_reference(m_wf, m_prop.m_ham, ref_onv, m_opts) {
    if (mpi::i_am_root())
        m_stats = mem_utils::make_unique<StatsFile<FciqmcStatsSpecifier>>("M7.stats");
    m_parallel_stats = mem_utils::make_unique<StatsFile<ParallelStatsSpecifier>>(
            "rank_"+std::to_string(mpi::irank())+".stats");
}

void Solver::execute(size_t niter) {
    for (size_t i = 0ul; i < niter; ++i) {
        reset();

        m_propagate_timer.unpause();
        loop_over_occupied_onvs();
        m_propagate_timer.pause();

        m_communicate_timer.unpause();
        m_wf.m_spawn.communicate();
        m_communicate_timer.pause();

        m_annihilate_timer.unpause();
        loop_over_spawned();
        m_annihilate_timer.pause();

        reduce();
        output_stats();
        ++m_icycle;
    }
}

void Solver::propagate_row(const size_t &irow) {
    if (m_wf.m_walkers.m_onv(irow).is_zero()) return;

    //const auto onv = m_wf.m_walkers.m_onv(irow);
    const auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
    if (consts::float_is_zero(weight)) return;

//    bool is_initiator = m_wf.m_walkers.m_flags.m_initiator(irow, 0, 0);
    bool is_deterministic = m_wf.m_walkers.m_flags.m_deterministic(irow);
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


    m_prop.off_diagonal(m_wf, irow);
    m_prop.diagonal(m_wf, irow);

    if (!is_deterministic && consts::float_is_zero(weight)) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "ALL WALKERS DIED: REMOVING DETERMINANT FROM LIST" << std::endl;
#endif
//                if (flag_initiator) m_ninitiator.m_delta--;
//                m_nocc_det.m_delta--;
//                m_data.remove(irow);
    }
}

void Solver::reset() {
    m_chk_nwalker_local = m_wf.m_nwalker(0, 0) + m_wf.m_delta_nwalker(0, 0);
    m_chk_ninitiator_local = m_wf.m_ninitiator(0, 0) + m_wf.m_delta_ninitiator(0, 0);
    m_wf.reset();
    m_reference.reset();
}

void Solver::reduce() {
    double chk_ratio;
    if (!consts::float_is_zero(m_wf.m_nwalker(0, 0))) {
        chk_ratio = m_chk_nwalker_local / m_wf.m_nwalker(0, 0);
        if (m_chk_nwalker_local > 0.0 && !consts::floats_nearly_equal(chk_ratio, 1.0))
            throw std::runtime_error("Unlogged walker population changes have occurred");
    }

    if (m_chk_ninitiator_local != m_wf.m_ninitiator(0, 0)) {
        throw std::runtime_error("Unlogged creations of initiator ONVs have occurred");
    }

    m_wf.reduce();
    m_reference.reduce();
    m_prop.update(m_icycle, m_wf);
}

void Solver::output_stats() {

    if (mpi::i_am_root()) {
        m_stats->m_icycle() = m_icycle;
        m_stats->m_tau() = m_prop.tau();
        m_stats->m_shift() = m_prop.m_shift;
        m_stats->m_nwalker() = m_wf.m_nwalker.reduced(0, 0);
        m_stats->m_delta_nwalker() = m_wf.m_delta_nwalker.reduced(0, 0);
        m_stats->m_ref_proj_energy_num() = m_reference.proj_energy_num();
        m_stats->m_ref_weight() = m_reference.weight();
        m_stats->m_ref_proj_energy() = m_reference.proj_energy();
        m_stats->m_l2_norm() = std::sqrt(m_wf.m_l2_norm_square.reduced(0, 0));
        m_stats->m_ninitiator() = m_wf.m_ninitiator.reduced(0, 0);
        m_stats->m_nocc_onv() = m_wf.m_nocc_onv.reduced(0, 0);
        m_stats->m_psingle() = m_prop.m_magnitude_logger.m_psingle;
        m_stats->flush();
    }


    m_parallel_stats->m_icycle() = m_icycle;
    m_parallel_stats->m_synchronization_wait_time() = m_synchronization_timer.total();
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
