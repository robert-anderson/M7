//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H

#include "src/data/PerforableMappedList.h"
#include <memory>

class WalkerList : public PerforableMappedList<Determinant> {
public:
    Field<Determinant> m_determinant;
    Field<defs::ham_t> m_weight;
    Field<defs::ham_comp_t> m_hdiag;
    Field<bool> m_reference_connection;
    Field<bool> m_initiator;
    Field<bool> m_deterministic;

    WalkerList(size_t nsite, size_t nrow);
};


#endif //M7_WALKERLIST_H
