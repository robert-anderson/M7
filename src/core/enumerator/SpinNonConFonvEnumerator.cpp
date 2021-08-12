//
// Created by rja on 04/11/2020.
//

#include "SpinNonConFonvEnumerator.h"

SpinNonConFonvEnumerator::SpinNonConFonvEnumerator(size_t nsite, size_t nelec) :
        Enumerator<field::FrmOnv>(), m_combs(2 * nsite, nelec), m_setinds(nelec, 0ul){}

bool SpinNonConFonvEnumerator::next_element(field::FrmOnv &result) {
    auto allfound = m_combs.next_element(m_setinds);
    if (!allfound) return false;
    result.zero();
    result = m_setinds;
    return true;
}
