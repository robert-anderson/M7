//
// Created by Robert John Anderson on 2020-04-03.
//

#include <src/core/thread/Atomic.h>
#include "Wavefunction.h"
#include "FciqmcCalculation.h"


Wavefunction::Wavefunction(FciqmcCalculation *fciqmc) :
    m_fciqmc(fciqmc), m_input(fciqmc->m_input), m_reference(m_fciqmc->m_reference),
    m_prop(fciqmc->m_prop), m_send(m_reference.nsite(), mpi::nrank()), m_recv(m_reference.nsite(), 1),
    m_data(m_reference.nsite(), m_input.nwalker_target * m_input.walker_factor_initial) {
    m_data.expand((size_t) (m_input.walker_factor_initial * m_input.nwalker_target));
    m_send.recv(&m_recv);
    m_send.expand((size_t) (m_input.buffer_factor_initial * m_input.nwalker_target));
    m_recv.expand((size_t) (m_input.buffer_factor_initial * m_input.nwalker_target));
    auto ref_weight = std::max(m_input.nwalker_initial, m_input.nadd_initiator);
    m_reference_row =
        m_data.add(m_reference, ref_weight, m_fciqmc->m_ham->get_energy(m_reference), true, true, false);
    m_square_norm = std::pow(std::abs(ref_weight), 2);
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

    m_data.synchronize();
    m_delta_square_norm = 0;
    m_delta_nw = 0;
    m_ref_proj_energy_num = 0;
    // capture the reference weight before death step is applied in the loop below
    m_reference_weight = *m_data.m_weight(m_reference_row);

#pragma omp parallel default(none)
    {
        defs::wf_comp_t delta_square_norm = 0;
        defs::wf_comp_t delta_nw = 0;
        defs::ham_t reference_energy_numerator = 0;
        int delta_ninitiator = 0;
#pragma omp for
        for (size_t irow = 0ul; irow < m_data.high_water_mark(0); ++irow) {
            if (m_data.row_empty(irow)) {
                continue;
            }
            auto weight = m_data.m_weight(irow);
            const auto det = m_data.m_determinant(irow);
            weight = m_prop->round(*weight);

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
                const auto contrib = *weight * m_prop->m_ham->get_element(m_reference, det);
                reference_energy_numerator += contrib;
            }
            auto hdiag = m_data.m_hdiag(irow);
            auto flag_deterministic = m_data.m_flags.m_deterministic(irow);

            //const_cast<DeterminantElement&>(det).print();
            m_prop->off_diagonal(det, weight, m_send, flag_deterministic, flag_initiator);
            m_prop->diagonal(hdiag, weight, delta_square_norm, delta_nw);
        }
        as_atomic(m_ninitiator) += delta_ninitiator;
        as_atomic(m_delta_square_norm) += delta_square_norm;
        as_atomic(m_delta_nw) += delta_nw;
        as_atomic(m_ref_proj_energy_num) += reference_energy_numerator;
    }
    assert(m_ninitiator >= 0);
}

void Wavefunction::communicate() {
    m_send.communicate();
}

void Wavefunction::annihilate_row(const size_t &irow_recv, defs::wf_comp_t &aborted_weight,
                                  defs::wf_comp_t &delta_square_norm, defs::wf_comp_t &delta_nw) {
    auto conn = m_fciqmc->m_scratch->conn->get(0);

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
        m_data.m_hdiag(irow_main) = m_prop->m_ham->get_energy(det);
        conn.zero();
        conn.connect(m_reference, det);
        if (conn.nexcit() < 3) {
            m_data.m_flags.m_reference_connection(irow_main) = true;
            assert(m_data.m_flags.m_reference_connection(irow_main));
        } else
            assert(!m_data.m_flags.m_reference_connection(irow_main));
    }
    auto weight = m_data.m_weight(irow_main);
    delta_square_norm += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
    delta_nw += std::abs(*weight + *delta_weight) - std::abs(*weight);
    weight += *delta_weight;
}

void Wavefunction::annihilate(){
    m_aborted_weight = 0;
#pragma omp parallel default(none)
    {
        Connection connection(m_reference);
        defs::wf_comp_t aborted_weight = 0;
        defs::wf_comp_t delta_square_norm = 0;
        defs::wf_comp_t delta_nw = 0;
#pragma omp for
        for (size_t irow_recv = 0ul; irow_recv < m_recv.high_water_mark(0); ++irow_recv) {
            annihilate_row(irow_recv, aborted_weight, delta_square_norm, delta_nw);
        }
        as_atomic(m_aborted_weight) += aborted_weight;
        as_atomic(m_delta_square_norm) += delta_square_norm;
        as_atomic(m_delta_nw) += delta_nw;
    }
    m_recv.zero();
    m_ninitiator = mpi::all_sum(m_ninitiator);
    m_delta_square_norm = mpi::all_sum(m_delta_square_norm);
    m_delta_nw = mpi::all_sum(m_delta_nw);
    m_square_norm = mpi::all_sum(m_square_norm);
    m_nw = mpi::all_sum(m_nw);
    m_noccupied_determinant = m_data.nfilled();
    m_noccupied_determinant = mpi::all_sum(m_noccupied_determinant);
    m_nw_growth_rate = (m_nw + m_delta_nw) / m_nw;
    m_square_norm += m_delta_square_norm;
    m_nw += m_delta_nw;
}

void Wavefunction::write_iter_stats(FciqmcStatsFile &stats_file) {
    //assert(m_data.m_flags.m_reference_connection(23));
    stats_file.m_ref_proj_energy_num() = m_ref_proj_energy_num;
    stats_file.m_ref_weight() = m_reference_weight;
    stats_file.m_ref_proj_energy() = m_ref_proj_energy_num / m_reference_weight;
    stats_file.m_nwalker() = m_nw;
    stats_file.m_aborted_weight() = m_aborted_weight;
    stats_file.m_ninitiator() = m_ninitiator;
    stats_file.m_noccupied_det() = m_noccupied_determinant;
}
