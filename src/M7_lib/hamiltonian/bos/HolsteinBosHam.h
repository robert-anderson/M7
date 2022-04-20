//
// Created by anderson on 12/9/21.
//

#ifndef M7_HOLSTEINBOSHAM_H
#define M7_HOLSTEINBOSHAM_H

#include "M7_lib/hamiltonian/bos/BosHam.h"

struct HolsteinBosHam : BosHam {
    defs::ham_comp_t m_omega;

    HolsteinBosHam(size_t nmode, const BosHilbertData& hd, defs::ham_comp_t omega):
        BosHam(BosBasisData{nmode}, hd), m_omega(omega){}

    defs::ham_t get_coeff_0011(size_t i, size_t j) const override {
        return i==j ? m_omega : 0;
    }

    defs::ham_t get_coeff_0022(size_t i, size_t j, size_t k, size_t l) const override {
        return 0;
    }

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override {
        defs::ham_t h = 0;
        for (size_t imode = 0ul; imode < m_bd.m_nmode; ++imode) h+= onv[imode];
        return h*m_omega;
    }

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0;
    }

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0;
    }
};


#endif //M7_HOLSTEINBOSHAM_H
