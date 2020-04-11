//
// Created by Robert John Anderson on 2020-04-03.
//

#include <src/core/thread/Atomic.h>
#include "Wavefunction.h"


Wavefunction::Wavefunction(const InputOptions &input, const std::unique_ptr<Propagator> &propagator,
                           const Determinant &reference) :
    m_input(input),
    m_send(reference.nsite(), mpi::nrank()), m_recv(reference.nsite(), 1),
    m_data(reference.nsite(), input.nwalker_target * input.walker_factor_initial),
    m_reference(reference) {
    m_data.expand((size_t) (input.walker_factor_initial * input.nwalker_target));
    m_send.recv(&m_recv);
    m_send.expand((size_t) (input.buffer_factor_initial * input.nwalker_target));
    m_recv.expand((size_t) (input.buffer_factor_initial * input.nwalker_target));
    auto ref_weight = std::max(m_input.nwalker_initial, m_input.nadd_initiator);
    m_reference_row =
        m_data.add(reference, ref_weight, propagator->m_ham->get_energy(reference), true, true, false);
    m_square_norm = std::pow(std::abs(ref_weight), 2);
}

void Wavefunction::propagate(std::unique_ptr<Propagator> &propagator) {
    /*
     * Loop over all rows in the m_data list and if the row is not empty:
     *      ascertain the initiator status of the determinant
     *      if not initiator and weight > initiator threshold, then grant initiator status
     *      if initiator and weight < initiator threshold, then revoke initiator status
     *
     *      perform all the off-diagonal propagation (fill send table)
     *      update local weight in the diagonal cloning/death step
     */

    m_data.synchronize();
    m_delta_square_norm = 0;
    m_ref_proj_energy_num = 0;
    // capture the reference weight before death step is applied in the loop below
    m_reference_weight = *m_data.m_weight(m_reference_row);

#pragma omp parallel default(none) shared(propagator)
    {
        defs::ham_comp_t delta_square_norm = 0;
        defs::ham_t reference_energy_numerator = 0;
        int delta_ninitiator = 0;
#pragma omp for
        for (size_t irow = 0ul; irow < m_data.high_water_mark(0); ++irow) {
            if (m_data.row_empty(irow)) {
                continue;
            }
            auto weight = m_data.m_weight(irow);
            const auto det = m_data.m_determinant(irow);
            weight = propagator->round(*weight);

            if (std::abs(*weight) == 0.0) {
                m_data.remove(det, irow);
                continue;
            }

            auto flag_initiator = m_data.m_flags.m_initiator(irow);
            if (!flag_initiator && std::abs(*weight) >= m_input.nadd_initiator) {
                // initiator status granted
                flag_initiator = true;
                delta_ninitiator++;
            } else if (flag_initiator && std::abs(*weight) < m_input.nadd_initiator) {
                // initiator status revoked
                flag_initiator = false;
                delta_ninitiator--;
            }
            if (m_data.m_flags.m_reference_connection(irow)) {
                const auto contrib = *weight * propagator->m_ham->get_element(m_reference, det);
                reference_energy_numerator += contrib;
            }
            auto hdiag = m_data.m_hdiag(irow);
            auto flag_deterministic = m_data.m_flags.m_deterministic(irow);

            //const_cast<DeterminantElement&>(det).print();
            propagator->off_diagonal(det, weight, m_send, flag_deterministic, flag_initiator);
            propagator->diagonal(hdiag, weight, delta_square_norm);
        }
        as_atomic(m_ninitiator) += delta_ninitiator;
        as_atomic(m_delta_square_norm) += delta_square_norm;
        as_atomic(m_ref_proj_energy_num) += reference_energy_numerator;
    }
    assert(m_ninitiator >= 0);
}

void Wavefunction::communicate() {
    m_send.communicate();
}

void Wavefunction::annihilate_row(const size_t &irow_recv, const std::unique_ptr<Propagator> &propagator,
                                  Connection &connection, defs::wf_comp_t &aborted_weight,
                                  defs::wf_comp_t &delta_square_norm) {
    auto det = m_recv.m_determinant(irow_recv);
    auto delta_weight = m_recv.m_weight(irow_recv);
    // zero magnitude weights should not have been communicated
    assert(!consts::float_is_zero(*delta_weight));
    size_t irow_main;

    auto mutex = m_data.key_mutex(det);
    irow_main = m_data.lookup(mutex, det);
    if (irow_main == ~0ul) {
        /*
         * the destination determinant is not currently occupied, so initiator rules
         * must be applied
         */
        if (!m_recv.m_flags.m_parent_initiator(irow_recv)) {
            aborted_weight += std::abs(*delta_weight);
            return;
        }
        irow_main = m_data.push(mutex, det);
        m_data.m_hdiag(irow_main) = propagator->m_ham->get_energy(det);
        connection.zero();
        connection.connect(m_reference, det);
        if (connection.nexcit() < 3) {
            m_data.m_flags.m_reference_connection(irow_main) = true;
            assert(m_data.m_flags.m_reference_connection(irow_main));
        }
        else assert(!m_data.m_flags.m_reference_connection(irow_main));
    }
    auto weight = m_data.m_weight(irow_main);
    delta_square_norm += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
    weight += *delta_weight;
}

void Wavefunction::annihilate(const std::unique_ptr<Propagator> &propagator) {
    m_aborted_weight = 0;
#pragma omp parallel default(none) shared(propagator)
    {
        Connection connection(m_reference);
        defs::wf_comp_t aborted_weight = 0;
        defs::wf_comp_t delta_square_norm = 0;
#pragma omp for
        for (size_t irow_recv = 0ul; irow_recv < m_recv.high_water_mark(0); ++irow_recv) {
            annihilate_row(irow_recv, propagator, connection, aborted_weight, delta_square_norm);
        }
        as_atomic(m_aborted_weight) += aborted_weight;
        as_atomic(m_delta_square_norm) += delta_square_norm;
    }
    m_recv.zero();
    m_ninitiator = mpi::all_sum(m_ninitiator);
    m_delta_square_norm = mpi::all_sum(m_delta_square_norm);
    m_square_norm = mpi::all_sum(m_square_norm);
    m_noccupied_determinant = m_data.nfilled();
    m_noccupied_determinant = mpi::all_sum(m_noccupied_determinant);
    //m_nw_growth_rate = std::sqrt((m_square_norm + m_delta_square_norm) / m_square_norm);
    m_square_norm += m_delta_square_norm;
    assert(m_data.m_flags.m_reference_connection(23));
}

defs::ham_comp_t Wavefunction::norm() const {
    return std::sqrt(m_square_norm);
}

void Wavefunction::write_iter_stats(FciqmcStatsFile &stats_file) {
    //assert(m_data.m_flags.m_reference_connection(23));
    stats_file.m_ref_proj_energy_num() = m_ref_proj_energy_num;
    stats_file.m_ref_weight() = m_reference_weight;
    stats_file.m_ref_proj_energy() = m_ref_proj_energy_num/m_reference_weight;
    stats_file.m_nwalker() = norm();
    stats_file.m_aborted_weight() = m_aborted_weight;
    stats_file.m_ninitiator() = m_ninitiator;
    stats_file.m_noccupied_det() = m_noccupied_determinant;
}