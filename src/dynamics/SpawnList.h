//
// Created by rja on 27/02/2020.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H


#include <cstddef>
#include "src/fermion/Determinant.h"
#include "src/data/Specification.h"
#include "src/data/List.h"

struct SpawnListFields {
    Specification m_spec;
    size_t idet;
    size_t iweight;
    size_t iflag_parent_initiator;

    SpawnListFields(size_t nsite) {
        idet = m_spec.add<Determinant>(nsite);
        iweight = m_spec.add<defs::ham_t>(1);
        iflag_parent_initiator = m_spec.add<bool>(1);
    }

};

class SpawnList : public List {
    SpawnListFields m_fields;
public:
    SpawnList(size_t nsite, size_t nrow, defs::data_t *data_external);

    SpawnList(size_t nsite, size_t nrow);

    Determinant get_determinant(const size_t &irow);

    NumericView<defs::ham_t> get_weight(const size_t &irow);

    NumericView<bool> get_flag_parent_initiator(const size_t &irow);

    size_t push(const Determinant &determinant, const defs::ham_t &weight,
                bool flag_parent_initiator) {
        auto irow = List::push();
        get_determinant(irow) = determinant;
        get_weight(irow) = weight;
        get_flag_parent_initiator(irow) = flag_parent_initiator;
    }
};


#endif //M7_SPAWNLIST_H
