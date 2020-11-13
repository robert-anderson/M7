//
// Created by Robert John Anderson on 2020-04-03.
//

#include <src/core/io/Logging.h>
#include "Wavefunction.h"
#include "FciqmcCalculation.h"

#if 0
Wavefunction::Wavefunction(FciqmcCalculation *fciqmc) :
        m_fciqmc(fciqmc), m_input(fciqmc->m_input),
        m_prop(fciqmc->m_prop),
        m_data("wavefunction walker list", fciqmc->m_reference.nsite(),
               m_input.nwalker_target * m_input.walker_factor_initial),
        m_send("wavefunction outgoing spawn list", fciqmc->m_reference.nsite(), mpi::nrank()),
        m_recv("wavefunction incoming spawn list", fciqmc->m_reference.nsite(), 1),
        m_reference(m_data, m_fciqmc->m_rank_allocator, fciqmc->m_reference, m_input.reference_redefinition_thresh){

    const auto nrow_walker = (size_t) (m_input.nwalker_target*m_input.walker_factor_initial);
    m_data.expand(nrow_walker);
    m_send.recv(&m_recv);
    const auto nrow_buffer = (size_t) (m_input.buffer_factor_initial * m_input.nwalker_target);
    m_send.expand(nrow_buffer);

    if (m_reference.is_mine()) {
        auto ref_weight = std::max(m_input.nwalker_initial, m_input.nadd_initiator);
        auto ref_energy = m_fciqmc->m_ham->get_energy(m_reference);
        logger::write("Reference energy: " + std::to_string(ref_energy));
        m_data.m_weight(m_reference.irow()) = ref_weight;
        m_data.m_hdiag(m_reference.irow()) = ref_energy;
        m_data.m_flags.m_reference_connection(m_reference.irow()) = true;
        m_data.m_flags.m_initiator(m_reference.irow()) = true;

        m_ninitiator = 1;
        m_square_norm = std::pow(std::abs(ref_weight), 2);
        m_nwalker = std::abs(ref_weight);
        m_nocc_det.local()++;
    } else {
        m_ninitiator = 0;
        m_nwalker = 0;
    }
}

void Wavefunction::update(const size_t &icycle) {
    /*
     * parallelization stats are written out on each node
     */
    parallel_stats_file()->m_walker_list_high_water_mark.write(m_data.high_water_mark(0));
    parallel_stats_file()->m_walker_list_high_water_mark_fraction.write(m_data.high_water_mark(0)/(double)m_data.nrow_per_segment());

    nrow_free = 0ul;

    m_nwalker.accumulate();
    ASSERT(m_nwalker.m_delta == 0.0)
    ASSERT(consts::floats_nearly_equal(m_nwalker.reduced()/nwalker_check(), 1.0))
    m_square_norm.accumulate();
    ASSERT(m_square_norm.m_delta == 0.0)
    //ASSERT(m_square_norm.reduced()==m_data.square_norm(0))
    m_ninitiator.accumulate();
    ASSERT(m_ninitiator.m_delta == 0.0)
    ASSERT(m_ninitiator.reduced()==ninitiator_check())
    m_nocc_det.accumulate();
    ASSERT(m_nocc_det.m_delta == 0)
    ASSERT(m_nocc_det.reduced()==nocc_check())

#ifndef NDEBUG
    size_t ninitiator_verify = m_data.verify_ninitiator(m_input.nadd_initiator);
    ninitiator_verify = mpi::all_sum(ninitiator_verify);
    ASSERT(m_ninitiator.reduced() == ninitiator_verify);
#endif

    m_reference.update();

    auto& shift_epoch = m_prop->variable_shift();
    auto& semistoch_epoch = m_prop->semi_stochastic();
    if (shift_epoch) {
        if (semistoch_epoch.update(icycle,
                m_input.do_semistochastic && (icycle > shift_epoch.start() + m_input.ncycle_init_detsub))){
            ASSERT(!m_detsub)
            m_detsub = std::unique_ptr<DeterministicSubspace>(new DeterministicSubspace(m_data));
            if (m_input.nadd_thresh_semistoch>0.0){
                m_detsub->build_from_nadd_thresh(m_input.nadd_thresh_semistoch, m_input.nadd_initiator, m_prop->m_ham.get());
            }
            else if (m_input.walker_fraction_semistoch<1.0){
                m_detsub->build_from_nw_fraction(m_input.walker_fraction_semistoch, m_nwalker.reduced(), m_prop->m_ham.get());
            }
            else {
                m_detsub->build_from_det_connections(m_reference, m_prop->m_ham.get(), 2);
            }
            //m_detsub->build_from_whole_walker_list(m_prop->m_ham.get());
        }
        if (!m_mk_sums) m_mk_sums = std::unique_ptr<KramersSectorOccupation>(new KramersSectorOccupation(m_reference));
    }

    if (consts::float_is_zero(m_nwalker.reduced())) throw (std::runtime_error("All walkers died."));
    ASSERT(consts::floats_nearly_equal(m_nwalker.reduced() / m_data.l1_norm(0), 1.0));

    // empty the list that receives new walkers
    m_recv.zero();
}


void Wavefunction::propagate() {
    /*
     * Loop over all rows in the m_data list and if the row is not empty:
     *      ascertain the initiator status of the determinant
     *      if not initiator and weight > initiator threshold, then grant initiator status
     *      if initiator and weight < initiator threshold, then revoke initiator status
     *
     *      perform all the off-diagonal propagation (fill send table)
     *      update local weight in the diagonal cloning/death step
     */
#ifdef VERBOSE_DEBUGGING
    m_data.print();
#endif
#ifndef NDEBUG
    auto chk_proj_energy = projected_energy_check(m_prop->m_ham.get(), m_reference);
#endif

    mpi::barrier(); m_propagation_timer.unpause();

    if (m_prop->semi_stochastic()) m_detsub->gather_and_project();

    for (size_t irow = 0ul; irow < m_data.high_water_mark(0); ++irow) {
        if (m_data.row_empty(irow)) {
            nrow_free++;
            continue;
        }
        auto weight = m_data.m_weight(irow);
        const auto det = m_data.m_determinant(irow);

#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "PARENT DETERMINANT" << std::endl;
        std::cout << consts::verb << "walker list row:  " << irow << std::endl;
        std::cout << consts::verb << "bitstring:        " << det.to_string() << std::endl;
        std::cout << consts::verb << "weight:           " << *weight << std::endl;
#endif
        //weight = m_prop->round(*weight);

        if (m_mk_sums) m_mk_sums->add(det, *weight);

        m_reference.log_candidate_weight(irow, std::abs(*weight));
        if (m_reference.in_redefinition_cycle()){
            /*
             * refresh reference connections
             */
            m_data.m_flags.m_reference_connection(irow) = m_reference.is_connected(det);
        }

        auto flag_initiator = m_data.m_flags.m_initiator(irow);

        if (consts::float_is_zero(*weight) && !m_data.m_flags.m_deterministic(irow)) {
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "ZERO WEIGHT: REMOVING FROM LIST" << std::endl;
            std::cout << consts::verb << "is initiator:     " << flag_initiator << std::endl;
            std::cout << consts::verb << "weight:           " << *weight << std::endl;
#endif
            if (flag_initiator) {
#ifdef VERBOSE_DEBUGGING
                std::cout << consts::verb << consts::chevs << "INITIATOR STATUS REVOKED: DETERMINANT REMOVED" << std::endl;
#endif
                m_ninitiator.m_delta--;
            }
            m_nocc_det.m_delta--;
            m_data.remove(irow);
            continue;
        }

        if (!flag_initiator && std::abs(*weight) >= m_input.nadd_initiator) {
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "INITIATOR STATUS GRANTED" << std::endl;
#endif
            flag_initiator = true;
            m_ninitiator.m_delta++;
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
        if (m_data.m_flags.m_reference_connection(irow)) {
            const auto contrib = *weight * m_prop->m_ham->get_element(m_reference, det);
            m_reference.proj_energy_num()+= contrib;
        }

        auto hdiag = m_data.m_hdiag(irow);
        auto flag_deterministic = m_data.m_flags.m_deterministic(irow);

        m_prop->off_diagonal(det, weight, m_send, flag_deterministic, flag_initiator);
        m_prop->diagonal(hdiag, weight, flag_deterministic,
                         m_square_norm.m_delta, m_nwalker.m_delta);

        if (!flag_deterministic && consts::float_is_zero(*weight)) {
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "ALL WALKERS DIED: REMOVING DETERMINANT FROM LIST" << std::endl;
#endif
            if (flag_initiator) m_ninitiator.m_delta--;
            m_nocc_det.m_delta--;
            m_data.remove(irow);
        }
    }
    mpi::barrier(); m_propagation_timer.pause();

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

#ifndef NDEBUG
    ASSERT(nrow_free-m_nocc_det.m_delta.local()==m_data.nzero_rows())

    auto chk_hwm = nrow_free+m_nocc_det.local();
    if (chk_hwm!=m_data.high_water_mark(0)) {
        m_data.print();
        auto chk_nrow_in_free_stack = m_data.nrow_in_free_stack();
        std::cout << "free rows in walker list " << chk_nrow_in_free_stack << std::endl;
    }
    ASSERT(chk_hwm==m_data.high_water_mark(0))
    ASSERT(chk_hwm-m_data.nzero_rows()==m_data.map_size())
    ASSERT(consts::floats_equal(chk_proj_energy, m_reference.proj_energy()))
#endif
}

void Wavefunction::communicate() {
    mpi::barrier(); m_communication_timer.unpause();

    const auto& hwm = m_send.high_water_mark();
    auto tot_sent = std::accumulate(hwm.begin(), hwm.end(), 0);
    parallel_stats_file()->m_nrow_sent.write(tot_sent);
    auto max_sent_iter = std::max_element(hwm.begin(), hwm.end());
    parallel_stats_file()->m_largest_nrow_sent.write(*max_sent_iter);
    parallel_stats_file()->m_largest_send_list_filled_fraction.write(*max_sent_iter/(double)m_send.nrow_per_segment());
    parallel_stats_file()->m_irank_largest_nrow_sent.write(std::distance(hwm.begin(), max_sent_iter));

    m_send.communicate();

    parallel_stats_file()->m_nrow_recv.write(m_recv.high_water_mark(0));
    parallel_stats_file()->m_nrow_recv.write(m_recv.high_water_mark(0)/(double)m_recv.nrow_per_segment());

    mpi::barrier(); m_communication_timer.pause();
}

void Wavefunction::annihilate_row(const size_t &irow_recv) {

    auto det = m_recv.m_determinant(irow_recv);
    // check that the received determinant has come to the right place
    ASSERT(m_fciqmc->m_rank_allocator.get_rank(det)==mpi::irank())
    auto delta_weight = m_recv.m_weight(irow_recv);
    // zero magnitude weights should not have been communicated
    ASSERT(!consts::float_is_zero(*delta_weight));
    size_t irow_main;

    irow_main = m_data.lookup_irow(det);

#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "bitstring:             " << det.to_string() << std::endl;
    std::cout << consts::verb << "delta weight:          " << *delta_weight << std::endl;
    std::cout << consts::verb << "found in walker list:  " << string_utils::yn(irow_main!=~0ul) << std::endl;
#endif

    if (irow_main == ~0ul) {
        /*
         * the destination determinant is not currently occupied, so initiator rules
         * must be applied
         */
        if (!m_recv.m_flags.m_parent_initiator(irow_recv)) {
            m_aborted_weight += std::abs(*delta_weight);
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs <<
                "ABORTED SPAWN TO UNOCCUPIED DETERMINANT: PARENT NON-INITIATOR" << std::endl;
#endif
            return;
        }
        irow_main = m_data.push(det);
        ASSERT(m_data.m_determinant(irow_main)==det)
        m_data.m_hdiag(irow_main) = m_prop->m_ham->get_energy(det);
        m_data.m_flags.m_reference_connection(irow_main) = m_reference.is_connected(det);
        m_data.m_flags.m_deterministic(irow_main) = false;
        m_data.m_flags.m_initiator(irow_main) = false;
        m_nocc_det.m_delta++;
    }

#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "row in walker list:    " << irow_main << std::endl;
#endif
    /*
     * if we have stochastically generated a connection between determinants in a deterministic
     * subspace, so we must reject this connection.
     */
    if (m_recv.m_flags.m_parent_deterministic(irow_recv) && m_data.m_flags.m_deterministic(irow_main)) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs <<
                  "ABORTED SPAWN: STOCHASTICALLY GENERATED A DETERMINISTIC CONNECTION" << std::endl;
#endif
        return;
    }
    auto weight = m_data.m_weight(irow_main);
    m_square_norm.m_delta += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
    m_nwalker.m_delta -= std::abs(*weight);
    weight += *delta_weight;
    m_nwalker.m_delta += std::abs(*weight);
}

void Wavefunction::annihilate() {
    mpi::barrier(); m_annihilation_timer.unpause();
    if (m_prop->semi_stochastic()) m_detsub->update_weights(m_prop->tau(), m_nwalker.m_delta);
    m_aborted_weight = 0;
    for (size_t irow_recv = 0ul; irow_recv < m_recv.high_water_mark(0); ++irow_recv) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "RECEIVED DETERMINANT" << std::endl;
        std::cout << consts::verb << "row in recv buffer:    " << irow_recv << std::endl;
#endif
        annihilate_row(irow_recv);
    }
    mpi::barrier(); m_annihilation_timer.pause();
#if 0
    m_data.clear_tombstones();
#endif
}

void Wavefunction::synchronize() {
    /*
     * This method handles the MPI communication involved in collating the Reducable members
     * of the Wavefunction class.
     * TODO: "communication syndicates" to avoid many blocking MPI reductions by combining in an array
     *
     * As with all iterative algorithms, it is easy to mix up which iteration corresponds
     * to the current state of a variable. In parallel algorithms this hazard is exacerbated.
     * Thus, all communications in this method are well-annotated to clarify this issue at every stage.
     *
     * In general however, variables which describe the "current state" of the wavefunction
     * (principal variables) refer to the pre-propagation state. This is because the current
     * state is only determinable by a deterministic loop over occupied determinants, which
     * only occurs once (loop_over_occupied_onvs() method) in the whole algorithm.
     *
     * Members of this class are conceptually divided according to the following scheme
     *
     * PRINCIPAL VARIABLES (not subject to propagation at synchronization)
     * ------------------------------------------------------------------------------------
     *
     *     ACCUMULATED PRINCIPAL VARIABLES (value persists between cycles and is updated
     *     by a delta variable in update())
     *     ---------------------------------------------------------------------------------
     *       - m_nocc_det
     *       - m_nwalker
     *       - m_square_norm
     *       - m_ninitiator
     *
     *     INSTANTANEOUS PRINCIPAL VARIABLES (value is reset to default value in update())
     *     ---------------------------------------------------------------------------------
     *       - m_reference_weight
     *       - m_ref_proj_energy_num
     *
     * DELTA VARIABLES (relates directly to the propagation just performed)
     * ------------------------------------------------------------------------------------
     *  - m_d_nwalker
     *  - m_d_square_norm
     *  - m_aborted_weight
     *
     *
     *  use the prefix "d_" to denote the delta variable containing the update to a principal
     *  variable. Delta variables do not persist between iterations and are always reset to
     *  their default value at update()
     *
     *  Any arithmetic expression involving a combination of principal and delta variables
     *  is a delta variable
     */

    /*
     *
     * m_nwalker contains the local and sum-reduced walker number (l1-norm of the wavefunction)
     * from the previous iteration (i.e. it has not been subject to the propagation just
     * performed)
     *
     * m_d_nwalker contains the thread-resolved change in m_nwalker for the propagation just performed,
     * this must be sum-reduced to set the local value
     */
    m_nwalker.reduce_delta();
    /*
     * local and reduced values of m_nwalker are not updated with those in m_d_nwalker at this stage,
     * because the statistics based on the principal variables e.g. the reference-projected energy
     * estimator are yet to be output.
     *
     * the delta due to the square_norm deserves the same treatment.
     */
    m_square_norm.reduce_delta();
    /*
     * m_aborted_weight is another thread-accumulated delta variable, but unlike the two norms, it
     * does not correspond to a principal variable. It is reduced here, consumed in write_iter_stats
     * and reset like any other delta variable in update()
     */
    m_aborted_weight.mpi_sum();

    /*
     * same for the change in occupied determinants
     */
    m_nocc_det.reduce_delta();

    m_ninitiator.reduce_delta();

    /*
     * the growth rate is the ratio of the new walker number (after update is called) to
     * the present value
     */
    m_nwalker_growth_rate = 1 + m_nwalker.m_delta.reduced() / m_nwalker.reduced();

    parallel_stats_file()->m_nrow_free_walker_list.write(nrow_free);
}

void Wavefunction::write_iter_stats(FciqmcStatsFile *stats_file) {
    if (!mpi::i_am_root()) return;
    stats_file->m_ref_proj_energy_num.write(m_reference.proj_energy_num().reduced());
    stats_file->m_ref_weight.write(m_reference.weight());
    stats_file->m_ref_proj_energy.write(m_reference.proj_energy());
    stats_file->m_nwalker.write(m_nwalker.reduced());
    stats_file->m_nw_growth_rate.write(m_nwalker_growth_rate);
    stats_file->m_nw_at_doubles.write(m_reference.nwalker_at_doubles().reduced()-std::abs(m_reference.weight()));
    stats_file->m_ref_candidate_weight.write(m_reference.candidate_weight().reduced());
    stats_file->m_aborted_weight.write(m_aborted_weight.reduced());
    stats_file->m_ninitiator.write(m_ninitiator.reduced());
    stats_file->m_noccupied_det.write(m_nocc_det.reduced());
    stats_file->m_prop_time.write(m_propagation_timer.lap());
    stats_file->m_comm_time.write(m_communication_timer.lap());
    stats_file->m_anni_time.write(m_annihilation_timer.lap());
}

ParallelizationStatsFile *Wavefunction::parallel_stats_file() {
    return m_fciqmc->m_parallel_stats_file.get();
}

#endif