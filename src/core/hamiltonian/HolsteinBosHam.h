//
// Created by anderson on 12/9/21.
//

#ifndef M7_HOLSTEINBOSHAM_H
#define M7_HOLSTEINBOSHAM_H

#include "BosHam.h"

struct HolsteinBosHam : BosHam {
    defs::ham_comp_t m_omega;
    HolsteinBosHam(size_t nmode, defs::ham_comp_t omega):
        BosHam(nmode, 0ul), m_omega(omega){}

    defs::ham_t get_coeff_0011(const size_t &i, const size_t &j) const override {
        return i==j ? m_omega : 0;
    }

    defs::ham_t get_coeff_0022(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override {
        return 0;
    }

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override {
        defs::ham_t h = 0;
        for (size_t imode = 0ul; imode < m_nmode; ++imode) {
            if (!onv[imode]) continue;
            defs::ham_comp_t occ = onv[imode];
            h += occ;
        }
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
