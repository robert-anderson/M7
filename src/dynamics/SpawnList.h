//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H

#include <memory>
#include "src/data/List.h"
#include "src/data/TableCommunicator.h"

struct SpawnListSpecification : public Specification {
    size_t idet;
    size_t iweight;
    size_t iflag_parent_initiator;
    SpawnListSpecification(size_t nsite) : Specification(){
        idet = add<Determinant>(nsite);
        iweight = add<defs::ham_t>(1);
        iflag_parent_initiator = add<bool>(1);
    }
};


class SpawnList : public TableCommunicator<List> {
    SpawnList(size_t nsite, size_t nrow_send, size_t nrow_recv):
    TableCommunicator<List>(SpawnListSpecification(nsite), nrow_send, nrow_recv);

    void create(const Determinant &det, const size_t &idst_rank, const defs::ham_t &weight, bool initiator){
        auto &list = m_send[idst_rank];
        size_t irow = list.push();
        list.view<Determinant>(irow, m_idet) = det;
        *list.view<defs::ham_t>(irow, m_iweight) = weight;
        *list.view<bool>(irow, m_iflag_parent_initiator) = initiator;
    }

    size_t nrecv() const{
        return m_recv.high_water_mark();
    }

    Determinant recv_det(const size_t &irow){
        return m_recv.view<Determinant>(irow, m_idet);
    }
    defs::ham_t recv_weight(const size_t &irow){
        return m_recv.view<defs::ham_t>(irow, m_iweight);
    }

};


#endif //M7_SPAWNLIST_H
