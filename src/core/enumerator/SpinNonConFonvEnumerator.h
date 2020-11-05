//
// Created by rja on 04/11/2020.
//

#ifndef M7_SPINNONCONFONVENUMERATOR_H
#define M7_SPINNONCONFONVENUMERATOR_H

#include "CombinationEnumerator.h"
#include "src/core/field/Views.h"

class SpinNonConFonvEnumerator : public Enumerator<views::FermionOnv> {
    CombinationEnumerator m_combs;
    defs::inds m_setinds;
public:
    SpinNonConFonvEnumerator(size_t nsite, size_t nelec);

    bool next_element(views::FermionOnv &result) override;
};


#endif //M7_SPINNONCONFONVENUMERATOR_H