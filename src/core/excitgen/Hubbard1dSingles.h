//
// Created by rja on 02/05/2021.
//

#ifndef M7_HUBBARD1DSINGLES_H
#define M7_HUBBARD1DSINGLES_H

#include "UniformSingles.h"

struct Hubbard1dSingles : public UniformSingles {
    const bool m_pbc;

    Hubbard1dSingles(const Hamiltonian& h, PRNG& prng):
        UniformSingles(h, prng), m_pbc(h.m_frm.m_model_attrs.is_hubbard_1d_pbc()){}

    bool draw(const size_t& exsig, const field::FrmOnv &src,
              CachedOrbs &orbs, defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const size_t& exsig, const field::FrmBosOnv &src_onv,
              CachedOrbs &orbs, defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override;

    size_t approx_nconn() const override;
};


#endif //M7_HUBBARD1DSINGLES_H
