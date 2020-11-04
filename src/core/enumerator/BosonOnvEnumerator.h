//
// Created by rja on 04/11/2020.
//

#ifndef M7_BOSONONVENUMERATOR_H
#define M7_BOSONONVENUMERATOR_H

#include "src/core/field/Views.h"
#include "ProductEnumerator.h"

class BosonOnvEnumerator : public Enumerator<views::BosonOnv> {
    ProductEnumerator m_prod;
    defs::inds m_setoccs;
public:
    BosonOnvEnumerator(size_t nmode, size_t occ_cutoff):
    Enumerator<views::BosonOnv>(), m_prod(nmode, occ_cutoff+1), m_setoccs(nmode){}
    bool next_element(views::BosonOnv& result) override {
        auto allfound = m_prod.next_element(m_setoccs);
        if(!allfound) return false;
        result = m_setoccs;
        return true;
    }
};

#endif //M7_BOSONONVENUMERATOR_H
