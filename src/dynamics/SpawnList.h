//
// Created by rja on 27/02/2020.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H


#include <cstddef>
#include "src/fermion/Determinant.h"
#include "src/data/Specification.h"
#include "src/data/List.h"

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

class SpawnList : public List {
public:
    using spec_T = SpawnListSpecification;
protected:
    spec_T m_spec;

public:
    SpawnList(const spec_T &spec, size_t nrow, defs::data_t *data_external);

    SpawnList(const spec_T &spec, size_t nrow);

    Determinant get_determinant(const size_t &irow);

    NumericView<defs::ham_t> get_weight(const size_t &irow);

    NumericView<bool> get_flag_parent_initiator(const size_t &irow);
};


#endif //M7_SPAWNLIST_H
