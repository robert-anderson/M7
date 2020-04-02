//
// Created by rja on 27/02/2020.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H

#if 0
#include <cstddef>
#include "src/fermion/Determinant.h"
#include "src/data/Specification.h"
#include "src/data/List.h"


class SpawnList : public List {
    Field<Determinant> m_determinant;
    Field<defs::ham_t> m_weight;
    Field<bool> m_parent_initiator;
public:
    SpawnList(size_t nsite, defs::data_t *data_external= nullptr);

    size_t push(const Determinant &determinant, const defs::ham_t &weight,
                bool flag_parent_initiator) {
        auto irow = List::push();
        m_determinant.get(irow) = determinant;
        m_weight.get(irow) = weight;
        m_parent_initiator.get(irow) = flag_parent_initiator;
    }
};


#endif //M7_SPAWNLIST_H
#endif //M7_SPAWNLIST_H
