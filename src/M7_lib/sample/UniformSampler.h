//
// Created by Robert J. Anderson on 04/08/2020.
//

#ifndef M7_UNIFORMSAMPLER_H
#define M7_UNIFORMSAMPLER_H

#if 0
#include "ExcitationGenerator.h"

class UniformSampler : public ExcitationGenerator {

public:
    UniformSampler(const FermionHamiltonian *h, PRNG &prng);

    bool draw_single(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     prob_t &prob, ham_t &helem, AntisymFermionOnvConnection &anticonn) override=0;

    bool draw_double(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, prob_t &prob, ham_t &helem,
                     AntisymFermionOnvConnection &anticonn) override =0;
};


#endif //M7_UNIFORMSAMPLER_H
#endif //M7_UNIFORMSAMPLER_H
