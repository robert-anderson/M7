//
// Created by rja on 26/08/2021.
//

#ifndef M7_FRMBOSUNIFORMKINETIC_H
#define M7_FRMBOSUNIFORMKINETIC_H

#include "ExcitGen.h"
#include "UniformSingles.h"

struct FrmBosUniformKinetic : public FrmBosExcitGen {
    UniformSingles m_singles;

    FrmBosUniformKinetic(const Hamiltonian &h, PRNG &prng);

    bool draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override;
};


#endif //M7_FRMBOSUNIFORMKINETIC_H
