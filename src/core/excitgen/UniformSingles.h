//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FermionExcitationGenerator.h"

class UniformSingles : public FermionExcitationGenerator {

public:
    UniformSingles(const Hamiltonian<>* ham, PRNG& prng);

    bool draw(const views::Onv<0> &src_onv, views::Onv<0> &dst_onv, const OccupiedOrbitals &occs,
                    const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
                    conn::Antisym<0> &anticonn) override;

};


#endif //M7_UNIFORMSINGLES_H
