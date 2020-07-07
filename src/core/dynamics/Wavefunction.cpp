//
// Created by Robert John Anderson on 2020-04-03.
//

#include <src/core/thread/Atomic.h>
#include <src/core/io/Logging.h>
#include "Wavefunction.h"
#include "FciqmcCalculation.h"


Wavefunction::Wavefunction(FciqmcCalculation *fciqmc) :
        m_fciqmc(fciqmc), m_input(fciqmc->m_input),
        m_prop(fciqmc->m_prop),
        m_data(fciqmc->m_reference.nsite(),
               m_input.nwalker_target * m_input.walker_factor_initial),
        m_send(fciqmc->m_reference.nsite(), mpi::nrank()),
        m_recv(fciqmc->m_reference.nsite(), 1),
        m_reference(m_data, m_fciqmc->m_rank_allocator, fciqmc->m_reference) {
    const auto nrow_walker = (size_t) (m_input.walker_factor_initial * m_input.nwalker_target);
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
    } else {
        m_ninitiator = 0;
        m_nwalker = 0;
    }
}


void Wavefunction::update(const size_t &icycle) {
    ASSERT(m_nwalker.m_delta.threads_zeroed())
    m_nwalker.accumulate();
    ASSERT(m_square_norm.m_delta.threads_zeroed())
    m_square_norm.accumulate();
    //ASSERT(m_square_norm.reduced()==m_data.square_norm(0))
    ASSERT(m_ninitiator.m_delta.threads_zeroed())
    m_ninitiator.accumulate();

#ifndef NDEBUG
    size_t ninitiator_verify = m_data.verify_ninitiator(m_input.nadd_initiator);
    ninitiator_verify = mpi::all_sum(ninitiator_verify);
    ASSERT(m_ninitiator.reduced() == ninitiator_verify);
#endif

    m_reference.update();

    auto tmp = m_prop->icycle_vary_shift();
    if (m_input.do_semistochastic && !m_in_semistochastic_epoch && tmp != ~0ul) {
        m_in_semistochastic_epoch = (icycle > tmp + m_input.niter_init_detsub);
        if (m_in_semistochastic_epoch) {
            std::cout << "Entering semistochastic epoch on MC cycle " << icycle << std::endl;
            m_detsub = std::unique_ptr<DeterministicSubspace>(new DeterministicSubspace(m_data));
            m_detsub->build_from_det_connections(m_reference, m_prop->m_ham.get());
            //m_detsub->build_from_whole_walker_list(m_prop->m_ham.get());
        }
    }

    if (consts::float_is_zero(m_nwalker.reduced())) throw (std::runtime_error("All walkers died."));
    ASSERT(consts::floats_nearly_equal(m_nwalker.reduced() / m_data.l1_norm(0), 1.0));

    // empty the list that receives new walkers
    m_recv.zero();
    // update space-recycling stacks in the sending PerforableMappedList
    m_data.synchronize();
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

    if (m_in_semistochastic_epoch) m_detsub->gather_and_project();

#pragma omp parallel for default(none) shared(stderr)
    for (size_t irow = 0ul; irow < m_data.high_water_mark(0); ++irow) {
        if (m_data.row_empty(irow)) {
            continue;
        }
        auto weight = m_data.m_weight(irow);
        const auto det = m_data.m_determinant(irow);

        //weight = m_prop->round(*weight);

        auto flag_initiator = m_data.m_flags.m_initiator(irow);

        if (consts::float_is_zero(*weight) && !m_data.m_flags.m_deterministic(irow)) {
            if (flag_initiator) m_ninitiator.m_delta.thread()--;
            m_data.remove(det, irow);
            continue;
        }

        if (!flag_initiator && std::abs(*weight) >= m_input.nadd_initiator) {
            // initiator status granted
            flag_initiator = true;
            m_ninitiator.m_delta.thread()++;
        }
        /*
        else if (flag_initiator && std::abs(*weight) < m_input.nadd_initiator) {
            // initiator status revoked
            // flag_initiator = false;
            // delta_ninitiator--;
        }
         */
        if (m_data.m_flags.m_reference_connection(irow)) {
            const auto contrib = *weight * m_prop->m_ham->get_element(m_reference, det);
            m_reference.proj_energy_num().thread() += contrib;
        }

        auto hdiag = m_data.m_hdiag(irow);
        auto flag_deterministic = m_data.m_flags.m_deterministic(irow);

        m_prop->off_diagonal(det, weight, m_send, flag_deterministic, flag_initiator);
        m_prop->diagonal(hdiag, weight, flag_deterministic,
                         m_square_norm.m_delta.thread(), m_nwalker.m_delta.thread());

        if (!flag_deterministic && consts::float_is_zero(*weight)) {
            if (flag_initiator) m_ninitiator.m_delta.thread()--;
#ifndef NDEBUG
            auto irow_removed = m_data.remove(det, irow);
            ASSERT(irow_removed == irow);
#else
            m_data.remove(det, irow);
#endif
        }
    }
}

void Wavefunction::communicate() {
    m_send.communicate();
}

void Wavefunction::annihilate_row(const size_t &irow_recv) {

    auto det = m_recv.m_determinant(irow_recv);
    auto delta_weight = m_recv.m_weight(irow_recv);
    // zero magnitude weights should not have been communicated
    ASSERT(!consts::float_is_zero(*delta_weight));
    size_t irow_main;

    auto mutex = m_data.key_mutex(det);
    irow_main = m_data.lookup(mutex, det);
    if (irow_main == ~0ul) {
        /*
         * the destination determinant is not currently occupied, so initiator rules
         * must be applied
         */
        if (!m_recv.m_flags.m_parent_initiator(irow_recv)) {
            m_aborted_weight.thread() += std::abs(*delta_weight);
            return;
        }
        irow_main = m_data.push(mutex, det);
        m_data.m_hdiag(irow_main) = m_prop->m_ham->get_energy(det);
        m_data.m_flags.m_reference_connection(irow_main) = m_reference.is_connected(det);
        m_data.m_flags.m_deterministic(irow_main) = false;
        m_data.m_flags.m_initiator(irow_main) = false;
    }
    /*
     * if we have stochastically generated a connection between determinants in a deterministic
     * subspace, so we must reject this connection.
     */
    if (m_recv.m_flags.m_parent_deterministic(irow_recv) && m_data.m_flags.m_deterministic(irow_main)) {
        return;
    }
    auto weight = m_data.m_weight(irow_main);
    m_square_norm.m_delta.thread() += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
    m_nwalker.m_delta.thread() -= std::abs(*weight);
    weight += *delta_weight;
    m_nwalker.m_delta.thread() += std::abs(*weight);
}

void Wavefunction::annihilate() {
    if (m_in_semistochastic_epoch) m_detsub->update_weights(m_prop->m_tau, m_nwalker.m_delta);
    m_aborted_weight = 0;
#pragma omp parallel for default(none)
    for (size_t irow_recv = 0ul; irow_recv < m_recv.high_water_mark(0); ++irow_recv) {
        annihilate_row(irow_recv);
    }
}

void Wavefunction::synchronize() {
    /*
     * This method handles the MPI communication involved in collating the Distributed and
     * Hybrid members of the Wavefunction class.
     * TODO: "communication syndicates" to avoid many blocking MPI reductions by combining in an array
     *
     * As with all iterative algorithms, it is easy to mix up which iteration corresponds
     * to the current state of a variable. In parallel and hybrid-parallel algorithms this
     * hazard is exacerbated. Thus, all communications in this method are well-annotated to
     * clarify this issue at every stage.
     *
     * In general however, variables which describe the "current state" of the wavefunction
     * (principal variables) refer to the pre-propagation state. This is because the current
     * state is only determinable by a deterministic loop over occupied determinants, which
     * only occurs once (propagate() method) in the whole algorithm.
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
    m_aborted_weight.put_thread_sum();
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
}

void Wavefunction::write_iter_stats(FciqmcStatsFile *stats_file) {
    if (!mpi::i_am_root()) return;
    stats_file->m_ref_proj_energy_num.write(m_reference.proj_energy_num().reduced());
    stats_file->m_ref_weight.write(m_reference.weight());
    stats_file->m_ref_proj_energy.write(m_reference.proj_energy());
    stats_file->m_nwalker.write(m_nwalker.reduced());
    stats_file->m_nw_growth_rate.write(m_nwalker_growth_rate);
    stats_file->m_aborted_weight.write(m_aborted_weight.reduced());
    stats_file->m_ninitiator.write(m_ninitiator.reduced());
    stats_file->m_noccupied_det.write(m_nocc_det.reduced());
}