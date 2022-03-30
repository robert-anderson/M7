//
// Created by rja on 26/08/2021.
//

#ifndef M7_LADDERHOPPINGPC_H
#define M7_LADDERHOPPINGPC_H

#include <M7_lib/sample/Aliaser.h>

#include "LadderHoppingUniform.h"

/**
 * pre-compute a one-per-node shared Aliaser akin to that of the HeatBathDoubles but this time for sampling the 1101 and
 * 1110 exsigs of the general fermion-boson Hamiltonian
 */
struct LadderHoppingPc : public LadderHoppingUniform {
    Aliaser m_pick_n_given_pq;

    LadderHoppingPc(const Hamiltonian &h, PRNG &prng);

    std::string description() const override;

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override;
};


#endif //M7_LADDERHOPPINGPC_H
