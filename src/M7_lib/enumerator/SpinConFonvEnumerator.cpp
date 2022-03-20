//
// Created by rja on 04/11/2020.
//

#include "SpinConFonvEnumerator.h"

SpinConFonvEnumerator::SpinConFonvEnumerator(size_t nsite, size_t nelec, int spin) :
        Enumerator<field::FrmOnv>(),
        m_alpha_comb(nsite, ci_utils::nalpha(nelec, spin)),
        m_beta_comb(nsite, ci_utils::nbeta(nelec, spin)),
        m_alpha_setinds(ci_utils::nalpha(nelec, spin)),
        m_beta_setinds(ci_utils::nbeta(nelec, spin))
{
    m_alpha_comb.next(m_alpha_setinds);
}

bool SpinConFonvEnumerator::next_element(field::FrmOnv &result) {
    bool inner_allfound = !m_beta_comb.next(m_beta_setinds);
    if (inner_allfound) {
        m_beta_comb.next(m_beta_setinds);
        if (!m_alpha_comb.next(m_alpha_setinds)) return false;
    }
    result.zero();
    result.set(m_alpha_setinds, m_beta_setinds);
    return true;
}
