//
// Created by Robert J. Anderson on 12/9/21.
//

#ifndef M7_NUMOPBOSHAM_H
#define M7_NUMOPBOSHAM_H

#include "M7_lib/hamiltonian/bos/BosHam.h"

/**
 * simply the bosonic number operator w * sum_i b_i+ b_i
 * where w is the m_weight scalar member
 */
struct NumOpBosHam : BosHam {
    const ham_comp_t m_weight;

    NumOpBosHam(const sys::bos::Basis& basis, ham_comp_t weight): BosHam(basis), m_weight(weight){}

    ham_t get_coeff_0011(uint_t a, uint_t i) const override {
        return a==i ? m_weight : 0;
    }

    ham_t get_coeff_0022(uint_t /*a*/, uint_t /*b*/, uint_t /*i*/, uint_t /*j*/) const override {
        return 0;
    }

    ham_t get_element_0000(const field::BosOnv& onv) const override {
        return onv.nboson() * m_weight;
    }

    ham_t get_element_0011(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const override {
        return 0;
    }

    ham_t get_element_0022(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const override {
        return 0;
    }
};


#endif //M7_NUMOPBOSHAM_H
