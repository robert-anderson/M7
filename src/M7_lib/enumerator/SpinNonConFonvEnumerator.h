//
// Created by rja on 04/11/2020.
//

#ifndef M7_SPINNONCONFONVENUMERATOR_H
#define M7_SPINNONCONFONVENUMERATOR_H

#include <M7_lib/field/Fields.h>

#include "CombinationEnumerator.h"

class SpinNonConFonvEnumerator : public Enumerator<field::FrmOnv> {
    CombinationEnumerator m_combs;
    defs::inds m_setinds;
public:
    SpinNonConFonvEnumerator(size_t nsite, size_t nelec);

    bool next_element(field::FrmOnv &result) override;
};


#endif //M7_SPINNONCONFONVENUMERATOR_H
