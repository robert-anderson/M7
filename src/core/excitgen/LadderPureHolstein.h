//
// Created by rja on 04/08/2021.
//

#ifndef M7_LADDERPUREHOLSTEIN_H
#define M7_LADDERPUREHOLSTEIN_H

#include "ExcitGen.h"
#include "LadderPureUniform.h"

struct LadderPureHolstein : public LadderPureUniform {
    LadderPureHolstein(const Hamiltonian &h, PRNG &prng) : LadderPureUniform(h, prng) {}

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    std::string description() const override;
};


#endif //M7_LADDERPUREHOLSTEIN_H
