//
// Created by rja on 30/08/2021.
//

#ifndef M7_LADDERPUREHOLSTEINZPM_H
#define M7_LADDERPUREHOLSTEINZPM_H

#include "LadderPureHolstein.h"

struct LadderPureHolsteinZpm : LadderPureHolstein {
    LadderPureHolsteinZpm(const Hamiltonian &h, PRNG &prng): LadderPureHolstein(h, prng){}

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
              conn::FrmBosOnv &conn) override;

    std::string description() const override;
};


#endif //M7_LADDERPUREHOLSTEINZPM_H
