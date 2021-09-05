//
// Created by rja on 26/08/2021.
//

#ifndef M7_LADDERHOPPINGUNIFORM_H
#define M7_LADDERHOPPINGUNIFORM_H

#include "ExcitGen.h"
#include "UniformSingles.h"

struct LadderHoppingUniform : public LadderExcitGen {
    UniformSingles m_singles;

    LadderHoppingUniform(const Hamiltonian &h, PRNG &prng);

    bool draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    std::string description() const override;

    size_t approx_nconn() const override;
};


#endif //M7_LADDERHOPPINGUNIFORM_H
