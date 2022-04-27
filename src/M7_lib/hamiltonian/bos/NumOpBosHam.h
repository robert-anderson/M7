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
    const defs::ham_comp_t m_weight;

    NumOpBosHam(const sys::bos::Basis& basis, defs::ham_comp_t weight):
            BosHam(basis), m_weight(weight){}

    defs::ham_t get_coeff_0011(size_t i, size_t j) const override {
        return i==j ? m_weight : 0;
    }

    defs::ham_t get_coeff_0022(size_t i, size_t j, size_t k, size_t l) const override {
        return 0;
    }

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override {
        return onv.nboson() * m_weight;
    }

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0;
    }

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0;
    }
};


#endif //M7_NUMOPBOSHAM_H
