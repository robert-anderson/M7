//
// Created by Robert John Anderson on 2020-02-04.
//

#include "Wavefunction.h"
#include "src/core/io/Logging.h"

#if 0

Wavefunction::Wavefunction(const InputOptions &input, const std::unique_ptr<Propagator> &propagator,
                           const Determinant &reference) :
        m_input(input),
        m_walker_list(reference.nspatorb(), (size_t) (input.nwalker_target * input.walker_factor_initial)),
        m_walker_communicator(reference.nspatorb(),
                              (size_t) (input.nwalker_target * input.buffer_factor_initial),
                              (size_t) (input.nwalker_target * input.buffer_factor_initial * mpi::nrank())),
        m_reference(reference) {
    auto irow = m_walker_list.push(reference);
    auto ref_weight = m_walker_list.get_weight(irow);
    ref_weight = std::max(m_input.nwalker_initial, m_input.nadd_initiator);
    m_walker_list.get_hdiag(irow) = propagator->m_ham->get_energy(reference);
    m_walker_list.get_flag_reference_connection(irow) = true;

    m_square_norm = std::pow(std::abs(*ref_weight), 2);
}

void Wavefunction::propagate(std::unique_ptr<Propagator> &propagator) {
    m_delta_square_norm = 0;
    m_reference_energy_numerator = 0;
    // capture the reference weight before death step is applied in the loop below
    m_reference_weight = 0;
    size_t irow = m_walker_list.lookup(m_reference);
    if (irow != ~0ul) m_reference_weight = *m_walker_list.get_weight(irow);

//#pragma omp parallel default(none) shared(propagator)
    {
        defs::ham_comp_t delta_square_norm = 0;
        defs::ham_t reference_energy_numerator = 0;
        int delta_ninitiator = 0;
//#pragma omp for
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(); ++irow) {
            if (m_walker_list.row_empty(irow)) continue;
            auto weight = m_walker_list.get_weight(irow);
            auto det = m_walker_list.get_determinant(irow);
            weight = propagator->round(*weight);
            if (std::abs(*weight==0.0)) m_walker_list.remove(det, irow);

            auto flag_initiator = m_walker_list.get_flag_initiator(irow);
            if (!*flag_initiator && std::abs(*weight) >= m_input.nadd_initiator) {
                // initiator status granted
                flag_initiator = true;
                delta_ninitiator++;
            } else if (*flag_initiator && std::abs(*weight) < m_input.nadd_initiator) {
                // initiator status revoked
                flag_initiator = false;
                delta_ninitiator--;
            }
            if (*m_walker_list.get_flag_reference_connection(irow)) {
                reference_energy_numerator += *weight *
                                              propagator->m_ham->get_element(m_reference, det);
            }
            auto hdiag = m_walker_list.get_hdiag(irow);
            auto flag_deterministic = m_walker_list.get_flag_deterministic(irow);
            auto &send = m_walker_communicator.m_send;
            propagator->off_diagonal(det, weight, flag_deterministic, flag_initiator, send);
            propagator->diagonal(hdiag, weight, delta_square_norm);
        }
#pragma omp atomic update
        m_ninitiator += delta_ninitiator;
#pragma omp atomic update
        m_delta_square_norm += delta_square_norm;
#pragma omp critical
        m_reference_energy_numerator += reference_energy_numerator;
    }
}

void Wavefunction::communicate() {
    m_walker_communicator.communicate();
}

void Wavefunction::annihilate(const std::unique_ptr<Propagator> &propagator) {
    m_aborted_weight = 0;
#pragma omp parallel default(none) shared(propagator)
    {
        defs::ham_comp_t aborted_weight = 0;
        defs::ham_comp_t delta_square_norm = 0;
        auto &recv = m_walker_communicator.m_recv;
#pragma omp for
        for (size_t irow_recv = 0ul; irow_recv < recv.high_water_mark(); ++irow_recv) {
            auto det = recv.get_determinant(irow_recv);
            auto delta_weight = recv.get_weight(irow_recv);
            {
                auto mutex = m_walker_list.key_mutex(det);
                size_t irow_main;
                irow_main = m_walker_list.lookup(mutex, det);
                if (irow_main == ~0ul) {
                    /*
                     * the destination determinant is not currently occupied, so initiator rules
                     * must be applied
                     */
                    if (!(bool) recv.get_flag_parent_initiator(irow_recv)) {
                        aborted_weight += std::abs(*delta_weight);
                        continue;
                    }
                    irow_main = m_walker_list.push(mutex, det);
                    m_walker_list.get_hdiag(irow_main) = propagator->m_ham->get_energy(det);
                    if (m_reference.nexcit(det) < 3) {
                        m_walker_list.get_flag_reference_connection(irow_main) = true;
                    }
                }
                auto weight = m_walker_list.get_weight(irow_main);
                delta_square_norm += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
                weight += delta_weight;
            }
        }
#pragma omp atomic update
        m_aborted_weight += aborted_weight;
#pragma omp atomic update
        m_delta_square_norm += delta_square_norm;
    }
    m_walker_communicator.m_recv.zero();
    m_ninitiator = mpi::all_sum(m_ninitiator);
    m_delta_square_norm = mpi::all_sum(m_delta_square_norm);
    m_square_norm = mpi::all_sum(m_square_norm);
    m_noccupied_determinant = m_walker_list.nfilled();
    m_noccupied_determinant = mpi::all_sum(m_noccupied_determinant);
    m_norm_growth_rate = std::sqrt((m_square_norm + m_delta_square_norm) / m_square_norm);
    m_square_norm += m_delta_square_norm;
}

defs::ham_comp_t Wavefunction::norm() const {
    return std::sqrt(m_square_norm);
}

void Wavefunction::write_iter_stats(FciqmcStatsFile &stats_file) {
    stats_file.m_reference_projected_energy_numerator->write(m_reference_energy_numerator);
    stats_file.m_reference_weight->write(m_reference_weight);
    stats_file.m_reference_energy->write(m_reference_energy_numerator / m_reference_weight);
    stats_file.m_wavefunction_l2_norm->write(norm());
    stats_file.m_aborted_weight->write(m_aborted_weight);
    stats_file.m_ninitiator->write(m_ninitiator);
    stats_file.m_noccupied_determinant->write(m_noccupied_determinant);
}
#endif