//
// Created by rja on 04/08/2021.
//

#ifndef M7_LADDERPUREHOLSTEIN_H
#define M7_LADDERPUREHOLSTEIN_H

#include "ExcitGen.h"
#include "LadderPureUniform.h"

struct LadderPureHolstein : public LadderPureUniform {
    LadderPureHolstein(const Hamiltonian &h, PRNG &prng, size_t exsig) : LadderPureUniform(h, prng, {exsig}) {
        REQUIRE_TRUE(h.m_ladder.get(), "holstein excit gen requires defined ladder hamiltonian");
        REQUIRE_TRUE(dynamic_cast<const HolsteinLadderHam*>(h.m_ladder.get()),
                     "holstein excit gen requires holstein hamiltonian");
        REQUIRE_EQ(h.m_ladder->m_bd.m_nsite, h.m_ladder->m_bd.m_nmode,
                        "holstein excit gen assumes one boson mode per fermion site");
    }

    std::string description() const override;
};

struct LadderHolsteinCre : public LadderPureHolstein {
    LadderHolsteinCre(const Hamiltonian &h, PRNG &prng) : LadderPureHolstein(h, prng, exsig_utils::ex_0010){}

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
                     defs::prob_t &prob, conn::FrmBosOnv &conn) override;
};

struct LadderHolsteinAnn : public LadderPureHolstein {
    LadderHolsteinAnn(const Hamiltonian &h, PRNG &prng) : LadderPureHolstein(h, prng, exsig_utils::ex_0001){}

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
                     defs::prob_t &prob, conn::FrmBosOnv &conn) override;
};

#endif //M7_LADDERPUREHOLSTEIN_H
