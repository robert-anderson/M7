//
// Created by rja on 04/11/2020.
//

#ifndef M7_SPINNONCONDETENUMERATOR_H
#define M7_SPINNONCONDETENUMERATOR_H

#include "CombinationEnumerator.h"
#include "src/core/field/Views.h"

class SpinNonConDetEnumerator : public Enumerator<views::Determinant> {
    CombinationEnumerator m_combs;
    defs::inds m_setinds;
public:
    SpinNonConDetEnumerator(size_t nsite, size_t nelec):
    Enumerator<views::Determinant>(), m_combs(2*nsite, nelec), m_setinds(nelec, 0ul){}

    bool next_element(views::Determinant &result) override {
        auto allfound = m_combs.next_element(m_setinds);
        if (!allfound) return false;
        result.zero();
        result.set(m_setinds);
        return true;
    }
};


#endif //M7_SPINNONCONDETENUMERATOR_H