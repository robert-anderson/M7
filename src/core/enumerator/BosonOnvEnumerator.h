//
// Created by rja on 04/11/2020.
//

#ifndef M7_BOSONONVENUMERATOR_H
#define M7_BOSONONVENUMERATOR_H

#include "ProductEnumerator.h"
#include "src/core/field/Fields.h"

class BosonOnvEnumerator : public Enumerator<field::BosOnv> {
    ProductEnumerator m_prod;
    defs::inds m_setoccs;
public:
    BosonOnvEnumerator(size_t nmode, size_t occ_cutoff);
    bool next_element(field::BosOnv& result) override;
};

#endif //M7_BOSONONVENUMERATOR_H
