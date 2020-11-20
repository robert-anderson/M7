//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FermionExcitationGenerator.h"

class UniformSingles : public FermionExcitationGenerator {

public:
    UniformSingles(const FermionHamiltonian* ham, PRNG& prng);

    bool draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv, const OccupiedOrbitals &occ,
                    const VacantOrbitals &vac, defs::prob_t &prob, defs::ham_t &helem,
                    conn::AsFermionOnv &anticonn) override;

};


#endif //M7_UNIFORMSINGLES_H
