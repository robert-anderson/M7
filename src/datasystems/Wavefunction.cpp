//
// Created by Robert John Anderson on 2020-02-04.
//

#include "Wavefunction.h"
/*
void Wavefunction::merge_recv_with_store(const Propagator &prop) {
    for (auto irow_recv{0ul}; irow_recv < m_recv_buffer->highwatermark()[0]; ++irow_recv) {
        auto det = m_recv_buffer->view<Determinant>(0, irow_recv, m_onv_buffer_entry);
        auto tmp = m_recv_buffer->view<defs::ham_t>(0, irow_recv, m_weight_buffer_entry)[0];
        auto irow_store = 0;//m_store->safe_push(0, det);
        auto delta_weight = m_recv_buffer->view<defs::ham_t>(0, irow_recv, m_weight_buffer_entry)[0];
        auto stored_weight = m_store->view<defs::ham_t>(0, irow_store, 0);
        if (consts::float_is_zero(*stored_weight)){
            if (m_hf_connection_buffer_flag->get(*m_recv_buffer, 0, irow_recv)) {
                m_hf_connection_store_flag->set(*m_store, 0, irow_store);
            }
        }
        m_delta_square_norm += std::pow(abs(*stored_weight + delta_weight), 2);
        m_delta_square_norm -= std::pow(abs(*stored_weight), 2);
        m_delta_component_norm+=consts::component_norm(*stored_weight + delta_weight);
        m_delta_component_norm-=consts::component_norm(*stored_weight);
        *stored_weight += delta_weight;
        auto hdiag = get_hdiag(irow_store);
        if (consts::float_is_zero(*hdiag)) *hdiag = prop.m_h.get_energy(det);
    }
    m_recv_buffer->zero();
}

void Wavefunction::death(const Propagator &prop) {
    for (auto irow{0ul}; irow < m_store->highwatermark()[0]; ++irow) {
        auto hdiag = *get_hdiag(irow);
        auto stored_weight = get_stored_weight(irow);
        auto delta_weight = -prop.tau * (hdiag - prop.m_shift) * (*stored_weight);
        m_delta_square_norm += std::pow(abs(*stored_weight + delta_weight), 2);
        m_delta_square_norm -= std::pow(abs(*stored_weight), 2);
        m_delta_component_norm+=consts::component_norm(*stored_weight + delta_weight);
        m_delta_component_norm-=consts::component_norm(*stored_weight);
        *stored_weight += delta_weight;
    }
}

*/