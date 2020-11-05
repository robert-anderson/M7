//
// Created by rja on 04/11/2020.
//

#include "SpinNonConDetEnumerator.h"

SpinNonConDetEnumerator::SpinNonConDetEnumerator(size_t nsite, size_t nelec) :
        Enumerator<views::Determinant>(), m_combs(2*nsite, nelec), m_setinds(nelec, 0ul){}

bool SpinNonConDetEnumerator::next_element(views::Determinant &result) {
    auto allfound = m_combs.next_element(m_setinds);
    if (!allfound) return false;
    result.zero();
    result.set(m_setinds);
    return true;
}
