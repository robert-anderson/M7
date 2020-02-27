//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H

#include "src/data/PerforableMappedList.h"
#include <memory>

struct WalkerListSpecification : public Specification {
    size_t idet;
    size_t iweight;
    size_t ihdiag;
    size_t iflag_reference_connection;
    size_t iflag_initiator;
    size_t iflag_deterministic;
    WalkerListSpecification(size_t nsite) : Specification(){
        idet = add<Determinant>(nsite);
        iweight = add<defs::ham_t>(1);
        ihdiag = add<defs::ham_comp_t>(1);
        iflag_reference_connection = add<bool>(1);
        iflag_initiator = add<bool>(1);
        iflag_deterministic = add<bool>(1);
    }
};

class WalkerList : public PerforableMappedList<Determinant> {

public:
    WalkerList(const size_t &nsite, const size_t &nrow) :
    PerforableMappedList<Determinant>(WalkerListSpecification(nsite), nrow, 0){}

    const WalkerListSpecification &spec() const{
        return static_cast<const WalkerListSpecification &>(m_spec);
    }
};


#endif //M7_WALKERLIST_H
