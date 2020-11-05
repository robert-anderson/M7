//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHSAMPLERS_H
#define M7_HEATBATHSAMPLERS_H

#if 0

#include <src/core/hamiltonian/FermionHamiltonian.h>
#include <src/core/sample/ExcitationGenerator.h>
#include "src/core/sample/Aliaser.h"

/*
 * precomputed sampler for doubles
 */

class HeatBathSamplers : public ExcitationGenerator {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathSamplers(const FermionHamiltonian *h, PRNG &prng);

    bool draw_single(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     defs::prob_t &prob, defs::ham_t &helem, AntisymConnection &anticonn) override;

    bool draw_double(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem,
                     AntisymConnection &anticonn) override;
};


#endif //M7_HEATBATHSAMPLERS_H
#endif //M7_HEATBATHSAMPLERS_H
