//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H

#include <memory>
#include "src/data/List.h"
#include "src/data/TableCommunicator.h"

class SpawnList {
public:
    std::unique_ptr<TableCommunicator<List>> m_conn = nullptr;
    Specification m_spec{};

    size_t m_idet;
    size_t m_iweight;
    size_t m_iflag_parent_initiator;

    SpawnList(const size_t &nspinorb, const size_t &nrow_send, const size_t &nrow_recv) {
        m_idet = m_spec.add<Determinant>(nspinorb);
        m_iweight = m_spec.add<defs::ham_t>(1);
        m_iflag_parent_initiator = m_spec.add<bool>(1);
        m_conn = std::make_unique<TableCommunicator<List>>(m_spec, nrow_send, nrow_recv);
    }

    bool communicate() {
        return m_conn->communicate();
    }

    void create(const Determinant &det, const size_t &idst_rank, const defs::ham_t &weight, bool initiator){
        auto &list = m_conn->m_send[idst_rank];
        size_t irow = list.push();
        list.view<Determinant>(irow, m_idet) = det;
        *list.view<defs::ham_t>(irow, m_iweight) = weight;
        *list.view<bool>(irow, m_iflag_parent_initiator) = initiator;
    }

    size_t nrecv() const{
        return m_conn->m_recv.high_water_mark();
    }

    Determinant recv_det(const size_t &irow){
        return m_conn->m_recv.view<Determinant>(irow, m_idet);
    }
    defs::ham_t recv_weight(const size_t &irow){
        return *m_conn->m_recv.view<defs::ham_t>(irow, m_iweight);
    }

};


#endif //M7_SPAWNLIST_H
