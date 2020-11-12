//
// Created by rja on 10/11/2020.
//

#include "Solver.h"

void Solver::propagate() {
    /*
     * Loop over all rows in the m_data list and if the row is not empty:
     *      ascertain the initiator status of the ONV
     *      if not initiator and weight > initiator threshold, then grant initiator status
     *      if initiator and weight < initiator threshold, then revoke initiator status
     *      perform all the off-diagonal propagation (fill send table)
     *      update local weight in the diagonal cloning/death step
     */
    size_t nfree_row_chk = 0ul;
    for (size_t irow = 0ul; irow < m_wf.m_walkers.m_hwm; ++irow) {

        if (m_wf.m_walkers.m_onv(irow).is_zero()) {
            // row is empty
            ++nfree_row_chk;
            continue;
        }

        const auto onv = m_wf.m_walkers.m_onv(irow);
        const auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);

        auto is_initiator = m_wf.m_walkers.m_flags.m_initiator(irow, 0, 0);
        auto is_deterministic = m_wf.m_walkers.m_flags.m_deterministic(irow);
        auto is_ref_connection = m_wf.m_walkers.m_flags.m_reference_connection(irow);

        if (consts::float_is_zero(weight) && !is_deterministic) {
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "ZERO WEIGHT: REMOVING FROM LIST" << std::endl;
                std::cout << consts::verb << "is initiator:     " << flag_initiator << std::endl;
                std::cout << consts::verb << "weight:           " << *weight << std::endl;
#endif
            if (is_initiator) {
#ifdef VERBOSE_DEBUGGING
                std::cout << consts::verb << consts::chevs << "INITIATOR STATUS REVOKED: DETERMINANT REMOVED" << std::endl;
#endif
                //m_ninitiator.m_delta--;
            }
            //m_nocc_det.m_delta--;
            //m_data.remove(irow);
            continue;
        }


        if (!is_initiator && std::abs(weight) >= m_opts.nadd_initiator) {
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "INITIATOR STATUS GRANTED" << std::endl;
#endif
            is_initiator = true;
            //m_ninitiator.m_delta++;
        }
        /*
        else if (flag_initiator && std::abs(*weight) < m_input.nadd_initiator) {
    #ifdef VERBOSE_DEBUGGING
                std::cout << consts::verb << consts::chevs << "INITIATOR STATUS REVOKED: WEIGHT FELL BELOW THRESHOLD MAGNITUDE" << std::endl;
    #endif
            // initiator status revoked
            // flag_initiator = false;
            // delta_ninitiator--;
        }
         */
        if (is_ref_connection) {
            const auto contrib = weight * m_prop.m_ham.get_element(m_reference, onv);
            m_reference.proj_energy_num() += contrib;
        }

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
    mpi::barrier();
    //m_propagation_timer.pause();

    /*
     * the numerator of the reference projected energy Rayleigh quotient is treated like a
     * delta variable, but is an instantaneous principal variable determined by the pre-
     * propagation walker distribution. This also enables the instantaneous Rayleigh quotient
     * to be computed. This is a pure convenience to aid a "quick look" at the energy estimator
     * in the stats file - the numerator and denominator should be independently analyzed to
     * obtain the estimate. All of this is handled in the Reference class, so defer to the
     * method implemented therein.
     */
    m_reference.synchronize();

#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << consts::chevs << "END OF PROPAGATION LOOP CHECKS" << std::endl;
std::cout << consts::verb << "free rows found in walker list:    " << nrow_free << std::endl;
std::cout << consts::verb << "free rows after propagation:       " << m_data.nzero_rows() << std::endl;
std::cout << consts::verb << "occupied determinants before loop: " << m_nocc_det.local() << std::endl;
std::cout << consts::verb << "delta in occupied determinants:    " << m_nocc_det.m_delta.local() << std::endl;
std::cout << consts::verb << "map size:                          " << m_data.map_size() << std::endl;
std::cout << consts::verb << "high water mark:                   " << m_data.high_water_mark(0) << std::endl;
#endif

#if 0
    #ifndef NDEBUG
        ASSERT(nrow_free - m_nocc_det.m_delta.local() == m_data.nzero_rows())

        auto chk_hwm = nrow_free + m_nocc_det.local();
        if (chk_hwm != m_data.high_water_mark(0)) {
            m_data.print();
            auto chk_nrow_in_free_stack = m_data.nrow_in_free_stack();
            std::cout << "free rows in walker list " << chk_nrow_in_free_stack << std::endl;
        }
        ASSERT(chk_hwm == m_data.high_water_mark(0))
        ASSERT(chk_hwm - m_data.nzero_rows() == m_data.map_size())
        ASSERT(consts::floats_equal(chk_proj_energy, m_reference.proj_energy()))
#endif
#endif


}

void Solver::annihilate_row(const size_t &irow_recv) {
    auto &recv = m_wf.m_spawn.recv();
    auto dst_onv = recv.m_dst_onv(irow_recv);
    std::cout << dst_onv.to_string() << std::endl;
    if(dst_onv.is_zero()) exit(0);
    ASSERT(!dst_onv.is_zero());
    // check that the received determinant has come to the right place
    ASSERT(m_wf.m_ra.get_rank(dst_onv) == mpi::irank())
    auto delta_weight = recv.m_delta_weight(irow_recv, 0, 0);
    // zero magnitude weights should not have been communicated
    ASSERT(!consts::float_is_zero(delta_weight));

    auto irow_walkers = *m_wf.m_walkers[dst_onv];

#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "bitstring:             " << det.to_string() << std::endl;
        std::cout << consts::verb << "delta weight:          " << *delta_weight << std::endl;
        std::cout << consts::verb << "found in walker list:  " << string_utils::yn(irow_main!=~0ul) << std::endl;
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

        irow_walkers = m_wf.add_walker(
                dst_onv,
                delta_weight,
                m_prop.m_ham.get_energy(dst_onv),
                m_reference.is_connected(dst_onv),
                delta_weight  >= m_opts.nadd_initiator);
        //m_nocc_det.m_delta++;
    }
    else {
        m_wf.m_walkers.m_weight(irow_walkers, 0, 0)+=delta_weight;
    }

#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "row in walker list:    " << irow_main << std::endl;
#endif
    /*
     * if we have stochastically generated a connection between determinants in a deterministic
     * subspace, so we must reject this connection.
     */
    if (recv.m_flags.m_src_deterministic(irow_recv) && m_wf.m_walkers.m_flags.m_deterministic(irow_walkers)) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs <<
                  "ABORTED SPAWN: STOCHASTICALLY GENERATED A DETERMINISTIC CONNECTION" << std::endl;
#endif
        return;
    }
    auto weight = m_wf.m_walkers.m_weight(irow_walkers, 0, 0);
    //m_square_norm.m_delta += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
    //m_nwalker.m_delta -= std::abs(*weight);
    weight += delta_weight;
    //m_nwalker.m_delta += std::abs(*weight);
}

void Solver::annihilate() {
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
