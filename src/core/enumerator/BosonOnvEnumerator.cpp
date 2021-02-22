//
// Created by rja on 04/11/2020.
//

#include "BosonOnvEnumerator.h"

BosonOnvEnumerator::BosonOnvEnumerator(size_t nmode, size_t occ_cutoff) :
        Enumerator<fields::BosonOnv>(), m_prod(nmode, occ_cutoff + 1), m_setoccs(nmode){}

bool BosonOnvEnumerator::next_element(fields::BosonOnv &result) {
    auto allfound = m_prod.next_element(m_setoccs);
    if(!allfound) return false;
    result = m_setoccs;
    return true;
}
