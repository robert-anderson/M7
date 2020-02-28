//
// Created by Robert John Anderson on 2020-02-04.
//

#include "Wavefunction.h"
#include "src/io/Logging.h"

Wavefunction::Wavefunction(const std::unique_ptr<Propagator> &propagator, const Determinant &reference,
                           double nwalker_initial, size_t nrow_walkers, size_t nrow_send, size_t nrow_recv) :
        m_walker_list(reference.nspatorb(), nrow_walkers),
        m_walker_communicator(reference.nspatorb(), nrow_send, nrow_recv),
        m_reference(reference){
    auto irow = m_walker_list.push(reference);
    m_walker_list.get_weight(irow) = nwalker_initial;
    m_walker_list.get_hdiag(irow) = propagator->m_ham->get_energy(reference);
    m_walker_list.get_flag_initiator(irow) = true;
    m_walker_list.get_flag_reference_connection(irow) = true;

    m_square_norm = std::pow(std::abs(nwalker_initial), 2);
}

void Wavefunction::propagate(const std::unique_ptr<Propagator> &propagator) {
    m_delta_square_norm = 0;
    m_reference_energy_numerator = 0;
    // capture the reference weight before death step is applied in the loop below
    m_reference_energy_denominator = 0;
    size_t irow = m_walker_list.lookup(m_reference);
    if (irow!=~0ul) m_reference_energy_denominator = *m_walker_list.get_weight(irow);

#pragma omp parallel default(none) shared(m_walker_list, propagator)
    {
        defs::ham_comp_t delta_square_norm = 0;
        defs::ham_t reference_energy_numerator = 0;
#pragma omp for
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(); ++irow) {
            if (m_walker_list.row_empty(irow)) continue;
            auto det = m_walker_list.get_determinant(irow);
            auto weight = m_walker_list.get_weight(irow);
            if (*m_walker_list.get_flag_reference_connection(irow)){
                reference_energy_numerator+=*weight*
                        propagator->m_ham->get_element(m_reference, det);
            }
            auto hdiag = m_walker_list.get_hdiag(irow);
            auto flag_deterministic = m_walker_list.get_flag_deterministic(irow);
            auto flag_initiator = m_walker_list.get_flag_initiator(irow);
            auto &send = m_walker_communicator.m_send;
            propagator->off_diagonal(det, weight, flag_deterministic, flag_initiator, send);
            propagator->diagonal(hdiag, weight, delta_square_norm);
        }
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
#pragma omp parallel default(none) shared(m_walker_communicator, m_delta_square_norm, propagator)
    {
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
                    irow_main = m_walker_list.push(mutex, det);
                    m_walker_list.get_hdiag(irow_main) = propagator->m_ham->get_energy(det);
                    if (m_reference.nexcit(det)<3) {
                        m_walker_list.get_flag_reference_connection(irow_main) = true;
                    }
                }
                auto weight = m_walker_list.get_weight(irow_main);
                delta_square_norm += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
                weight += delta_weight;
            }
        }
#pragma omp atomic update
        m_delta_square_norm += delta_square_norm;
    }
    m_walker_communicator.m_recv.zero();
    m_delta_square_norm = mpi::all_sum(m_delta_square_norm);
    m_square_norm = mpi::all_sum(m_square_norm);
    m_norm_growth_rate = std::sqrt((m_square_norm+m_delta_square_norm) / m_square_norm);
    m_square_norm += m_delta_square_norm;
    logger::write("norm: "+std::to_string(std::sqrt(m_square_norm)));
}

defs::ham_comp_t Wavefunction::norm() const {
    return std::sqrt(m_square_norm);
}
