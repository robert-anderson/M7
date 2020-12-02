//
// Created by rja on 04/11/2020.
//

#ifndef M7_SPINCONFONVENUMERATOR_H
#define M7_SPINCONFONVENUMERATOR_H


#include "CombinationEnumerator.h"
#include "src/core/field/Views.h"

class SpinConFonvEnumerator : public Enumerator<views::Det> {
    CombinationEnumerator m_alpha_comb;
    CombinationEnumerator m_beta_comb;
    defs::inds m_alpha_setinds;
    defs::inds m_beta_setinds;
public:
    SpinConFonvEnumerator(size_t nsite, size_t nelec, int spin);

    bool next_element(views::Det &result) override;
};


#endif //M7_SPINCONFONVENUMERATOR_H
