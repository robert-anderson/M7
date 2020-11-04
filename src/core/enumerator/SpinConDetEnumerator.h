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
    SpinConDetEnumerator(size_t nsite, size_t nelec, int spin) :
            Enumerator<views::Determinant>(),
            m_alpha_comb(nsite, ci_utils::nalpha(nelec, spin)),
            m_beta_comb(nsite, ci_utils::nbeta(nelec, spin)),
            m_alpha_setinds(ci_utils::nalpha(nelec, spin)),
            m_beta_setinds(ci_utils::nbeta(nelec, spin))
    {
        m_alpha_comb.next(m_alpha_setinds);
    }

    bool next_element(views::Determinant &result) override {
        auto b_allfound = !m_beta_comb.next(m_beta_setinds);
        result.zero();
        result.set(m_alpha_setinds, m_beta_setinds);
        if (b_allfound) {
            auto a_allfound = !m_alpha_comb.next(m_alpha_setinds);
            if (a_allfound) return false;
        }
        return true;
    }
};


#endif //M7_SPINCONDETENUMERATOR_H
