//
// Created by rja on 04/11/2020.
//

#ifndef M7_SPINCONDETENUMERATOR_H
#define M7_SPINCONDETENUMERATOR_H


#include "CombinationEnumerator.h"
#include "src/core/field/Views.h"

class SpinConDetEnumerator : public Enumerator<views::Determinant> {
    CombinationEnumerator m_alpha_comb;
    CombinationEnumerator m_beta_comb;
    defs::inds m_alpha_setinds;
    defs::inds m_beta_setinds;
public:
    SpinConDetEnumerator(size_t nsite, size_t nelec, int spin);

    bool next_element(views::Determinant &result) override;
};


#endif //M7_SPINCONDETENUMERATOR_H
