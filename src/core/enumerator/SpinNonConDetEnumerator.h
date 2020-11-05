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
    SpinNonConDetEnumerator(size_t nsite, size_t nelec);

    bool next_element(views::Determinant &result) override;
};


#endif //M7_SPINNONCONDETENUMERATOR_H