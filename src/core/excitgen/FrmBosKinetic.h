//
// Created by rja on 26/08/2021.
//

#ifndef M7_FRMBOSKINETIC_H
#define M7_FRMBOSKINETIC_H

#include "src/core/sample/Aliaser.h"
#include "ExcitGen.h"
#include "UniformSingles.h"

/**
 * set up a one-per-node shared Aliaser akin to that of the HeatBathDoubles but this time for sampling the 1101 and 1101
 * exsigs of the general fermion-boson Hamiltonian
 */
struct FrmBosKinetic : public FrmBosExcitGen {
    Aliaser m_pick_n_given_pq;
    UniformSingles m_singles;

    FrmBosKinetic(const Hamiltonian &h, PRNG &prng);

    bool draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override;
};


#endif //M7_FRMBOSKINETIC_H
