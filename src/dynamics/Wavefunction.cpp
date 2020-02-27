//
// Created by Robert John Anderson on 2020-02-04.
//

#include "Wavefunction.h"



Wavefunction::Wavefunction(const Determinant &ref, size_t nwalker_initial, size_t nrow_walkers, size_t nrow_send,
                           size_t nrow_recv) :
        m_walker_list(ref.nspatorb(), nrow_walkers),
        m_walker_communicator(ref.nspatorb(), nrow_send, nrow_recv){
    auto irow = m_walker_list.push(ref);
    m_walker_list.get_weight(irow) = nwalker_initial;
    m_square_norm = std::pow(nwalker_initial, 2);
}

void Wavefunction::propagate(const std::unique_ptr<Propagator> &propagator) {
    m_delta_square_norm = 0;
#pragma omp parallel default(none) shared(m_walker_list, propagator)
    {
        defs::ham_comp_t delta_square_norm = 0;
#pragma omp for
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(); ++irow) {
            if (m_walker_list.row_empty(irow)) continue;
            auto det = m_walker_list.get_determinant(irow);
            auto weight = m_walker_list.get_weight(irow);
            auto hdiag = m_walker_list.get_hdiag(irow);
            auto flag_deterministic = m_walker_list.get_flag_deterministic(irow);
            auto flag_initiator = m_walker_list.get_flag_initiator(irow);
            auto &send = m_walker_communicator.m_send;
            propagator->off_diagonal(det, weight, flag_deterministic, flag_initiator, send);
            propagator->diagonal(hdiag, weight, delta_square_norm);
        }
#pragma omp atomic update
        m_delta_square_norm += delta_square_norm;
    }
}

void Wavefunction::communicate() {
    m_walker_communicator.communicate();
}

void Wavefunction::annihilate() {
#pragma omp parallel default(none) shared(m_walker_communicator, m_delta_square_norm)
    {
        defs::ham_comp_t delta_square_norm = 0;
        auto &recv = m_walker_communicator.m_recv;
        for (size_t irow_recv = 0ul; irow_recv < recv.high_water_mark(); ++irow_recv) {
            auto det = recv.get_determinant(irow_recv);
            auto delta_weight = recv.get_weight(irow_recv);
            {
                auto mutex = m_walker_list.key_mutex(det);
                size_t irow_main = m_walker_list.push(mutex, det);
                auto weight = m_walker_list.get_weight(irow_main);
                delta_square_norm+=std::pow(std::abs(*weight+*delta_weight), 2)-std::pow(std::abs(*weight), 2);
                weight+=delta_weight;
            }
        }
#pragma omp atomic update
        m_delta_square_norm += delta_square_norm;
    }
    m_walker_communicator.m_recv.zero();
    m_delta_square_norm = mpi::all_sum(m_delta_square_norm);
    m_square_norm = mpi::all_sum(m_square_norm);
    m_norm_growth_rate = 1.0+std::sqrt(m_delta_square_norm/m_square_norm);
    m_square_norm+=m_delta_square_norm;
}
