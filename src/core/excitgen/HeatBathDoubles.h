//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHDOUBLES_H
#define M7_HEATBATHDOUBLES_H

#include <src/core/hamiltonian/FermionHamiltonian.h>
#include <src/core/excitgen/FermionExcitationGenerator.h>
#include "src/core/sample/Aliaser.h"
#include "src/core/field/Views.h"

/*
 * precomputed sampler for doubles
 */

class HeatBathDoubles : public FermionExcitationGenerator {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const FermionHamiltonian *h, PRNG &prng);

    bool draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     defs::prob_t &prob, defs::ham_t &helem, conn::AsFermionOnv &anticonn) override;

};

#endif //M7_HEATBATHDOUBLES_H
