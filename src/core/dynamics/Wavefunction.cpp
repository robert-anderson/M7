#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by Robert John Anderson on 2020-04-03.
//

#include <src/core/thread/Atomic.h>
#include <src/core/io/Logging.h>
#include "Wavefunction.h"
#include "FciqmcCalculation.h"


Wavefunction::Wavefunction(FciqmcCalculation *fciqmc) :
        m_fciqmc(fciqmc), m_input(fciqmc->m_input), m_reference(m_fciqmc->m_reference),
        m_prop(fciqmc->m_prop), m_send(m_reference.nsite(), mpi::nrank()), m_recv(m_reference.nsite(), 1),
        m_data(m_reference.nsite(), m_input.nwalker_target * m_input.walker_factor_initial) {
    const auto nrow_walker = (size_t) (m_input.walker_factor_initial * m_input.nwalker_target);
    m_data.expand(nrow_walker);
    m_send.recv(&m_recv);
    const auto nrow_buffer = (size_t) (m_input.buffer_factor_initial * m_input.nwalker_target);
    m_send.expand(nrow_buffer);
    auto ref_weight = std::max(m_input.nwalker_initial, m_input.nadd_initiator);
    auto ref_energy = m_fciqmc->m_ham->get_energy(m_reference);
    logger::write("Reference energy: " + std::to_string(ref_energy));

    m_irank_reference = fciqmc->m_rank_allocator.get_rank(m_reference);
    m_irank_reference.mpi_bcast(0);
    if (mpi::i_am(m_irank_reference.reduced())) {
        m_reference_row =
                m_data.add(m_reference, ref_weight, ref_energy, true, true, false);
        m_ninitiator.local() = 1;
        m_square_norm = std::pow(std::abs(ref_weight), 2);
        m_nw = std::abs(ref_weight);
    } else {
        m_reference_row = ~0ul;
        m_ninitiator.local() = 0;
        m_nw = 0;
    }
    m_nw.mpi_sum();
    m_square_norm.mpi_sum();
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
    m_ref_proj_energy_num = 0;

#pragma omp parallel for
    for (size_t irow = 0ul; irow < m_data.high_water_mark(0); ++irow) {
        if (m_data.row_empty(irow)) {
            continue;
        }
        auto weight = m_data.m_weight(irow);
        const auto det = m_data.m_determinant(irow);
        weight = m_prop->round(*weight);
        auto flag_initiator = m_data.m_flags.m_initiator(irow);

        if (consts::float_is_zero(*weight)) {
            if (flag_initiator) m_ninitiator.thread()--;
            m_data.remove(det, irow);
            continue;
        }

        if (!flag_initiator && std::abs(*weight) >= m_input.nadd_initiator) {
            // initiator status granted
            flag_initiator = true;
            m_ninitiator.thread()++;
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
            m_ref_proj_energy_num.thread() += contrib;
        }
        auto hdiag = m_data.m_hdiag(irow);
        auto flag_deterministic = m_data.m_flags.m_deterministic(irow);

        m_prop->off_diagonal(det, weight, m_send, flag_deterministic, flag_initiator);
        m_prop->diagonal(hdiag, weight, m_delta_square_norm.thread(), m_delta_nw.thread());

        if (consts::float_is_zero(*weight)) {
            if (flag_initiator) m_ninitiator.thread()--;
            auto irow_removed = m_data.remove(det, irow);
            ASSERT(irow_removed == irow);
        }
    }
    m_ninitiator.add_thread_sum();
    m_ref_proj_energy_num.put_thread_sum();
}

void Wavefunction::communicate() {
    m_send.communicate();
}

void Wavefunction::annihilate_row(const size_t &irow_recv) {

    auto conn = m_fciqmc->m_scratch->conn->get();
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
        conn.zero();
        conn.connect(m_reference, det);
        m_data.m_flags.m_reference_connection(irow_main) = conn.nexcit() < 3;
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
    m_delta_square_norm.thread() += std::pow(std::abs(*weight + *delta_weight), 2) - std::pow(std::abs(*weight), 2);
    m_delta_nw.thread() += std::abs(*weight + *delta_weight) - std::abs(*weight);
    weight += *delta_weight;
}

void Wavefunction::annihilate() {
    m_aborted_weight = 0;
#pragma omp parallel for
    for (size_t irow_recv = 0ul; irow_recv < m_recv.high_water_mark(0); ++irow_recv) {
        annihilate_row(irow_recv);
    }
    m_aborted_weight.put_thread_sum();
    m_delta_square_norm.put_thread_sum();
    m_delta_nw.put_thread_sum();
    m_recv.zero();
    m_ninitiator.mpi_sum();
    ASSERT(m_ninitiator.reduced() >= 0)
#ifndef NDEBUG
    if (mpi::nrank() == 1) {
        size_t ninitiator_verify = m_data.verify_ninitiator(m_input.nadd_initiator);
        DBVAR(ninitiator_verify);
        DBVAR(m_ninitiator.reduced());
        ASSERT(m_ninitiator.reduced() == ninitiator_verify);
    }
#endif

    m_delta_square_norm.mpi_sum();
    m_delta_nw.mpi_sum();
    m_square_norm.mpi_sum();
    m_noccupied_determinant.mpi_sum();
    m_nw_growth_rate = (m_nw.reduced() + m_delta_nw.reduced()) / m_nw.reduced();
    m_square_norm.local() += m_delta_square_norm.local();
    m_nw.local()+=m_delta_nw.local();
    m_nw.mpi_sum();
    if (consts::float_is_zero(m_nw.reduced())) throw (std::runtime_error("All walkers died."));
    m_ref_proj_energy_num.mpi_sum();
    m_aborted_weight.mpi_sum();

    if (mpi::i_am(m_irank_reference.reduced())) {
        m_reference_weight = *m_data.m_weight(m_reference_row);
    }
    m_reference_weight.mpi_bcast(m_irank_reference.reduced());
}

void Wavefunction::write_iter_stats(FciqmcStatsFile *stats_file) {
    if (!mpi::i_am_root()) return;
    stats_file->m_ref_proj_energy_num.write(m_ref_proj_energy_num.reduced());
    stats_file->m_ref_weight.write(m_reference_weight.reduced());
    stats_file->m_ref_proj_energy.write(m_ref_proj_energy_num.reduced() / m_reference_weight.reduced());
    stats_file->m_nwalker.write(m_nw.reduced());
    stats_file->m_aborted_weight.write(m_aborted_weight.reduced());
    stats_file->m_ninitiator.write(m_ninitiator.reduced());
    stats_file->m_noccupied_det.write(m_noccupied_determinant.reduced());
}

#pragma clang diagnostic pop