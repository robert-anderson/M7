//
// Created by rja on 02/05/2021.
//

#ifndef M7_HUBBARDUNIFORM_H
#define M7_HUBBARDUNIFORM_H

#include <hamiltonian/HubbardFrmHam.h>

#include "UniformSingles.h"

struct HubbardUniform : public UniformSingles {
    using UniformSingles::draw;

    const HubbardFrmHam* h_cast() {
        return dynamic_cast<const HubbardFrmHam*>(m_h.m_frm.get());
    }

    HubbardUniform(const Hamiltonian& h, PRNG& prng): UniformSingles(h, prng) {
        REQUIRE_TRUE(h_cast(), "fermion hamiltonian is not of hubbard hamiltonian type");
    }

    bool draw_frm(const size_t& exsig, const field::FrmOnv &src,
                  CachedOrbs &orbs, defs::prob_t &prob, conn::FrmOnv &conn) override;

    bool draw_frmbos(const size_t& exsig, const field::FrmBosOnv &src_onv,
              CachedOrbs &orbs, defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    size_t approx_nconn() const override;

    std::string description() const override;

};


#endif //M7_HUBBARDUNIFORM_H
