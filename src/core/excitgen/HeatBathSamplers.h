//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHSAMPLERS_H
#define M7_HEATBATHSAMPLERS_H

#include <src/core/hamiltonian/FermionHamiltonian.h>
#include <src/core/sample/ExcitationGenerator.h>
#include "src/core/sample/Aliaser.h"
#include "src/core/field/Views.h"

/*
 * precomputed sampler for doubles
 */

class HeatBathSamplers : public ExcitationGenerator {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathSamplers(const FermionHamiltonian *h, PRNG &prng);

    bool draw_single(const views::FermionOnv &src_fonv, views::Onv &dst_fonv,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     defs::prob_t &prob, defs::ham_t &helem, conn::AsFermionOnv &anticonn) override;

    bool draw_double(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem,
                     conn::AsFermionOnv &anticonn) override;
};

#endif //M7_HEATBATHSAMPLERS_H
