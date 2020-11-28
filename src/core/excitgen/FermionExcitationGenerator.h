//
// Created by RJA on 20/11/2020.
//

#ifndef M7_FERMIONEXCITATIONGENERATOR_H
#define M7_FERMIONEXCITATIONGENERATOR_H

#include "ExcitationGenerator.h"

class FermionExcitationGenerator : public ExcitationGenerator {
protected:
    const size_t m_nexcit;
    const bool m_spin_conserving;

public:
    FermionExcitationGenerator(const Hamiltonian<> *h, PRNG &prng, size_t nexcit);

    bool draw(const views::Onv<0> &src_onv, views::Onv<0> &dst_onv,
              const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
              defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<0> &anticonn) override;

    bool draw(const views::Onv<1> &src_onv, views::Onv<1> &dst_onv,
              const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
              defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<1> &anticonn) override;
};


#endif //M7_FERMIONEXCITATIONGENERATOR_H
