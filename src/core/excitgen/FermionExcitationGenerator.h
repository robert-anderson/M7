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
    FermionExcitationGenerator(const Hamiltonian *h, PRNG &prng, size_t nexcit);

    bool draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv, const OccupiedOrbitals &occ,
                    const VacantOrbitals &vac, defs::prob_t &prob, defs::ham_t &helem,
                    conn::AsFermionOnv &anticonn) override;

    bool draw(const views::FermiBosOnv &src_fonv, views::FermiBosOnv &dst_fonv, const OccupiedOrbitals &occ,
                    const VacantOrbitals &vac, defs::prob_t &prob, defs::ham_t &helem,
                    conn::AsFermiBosOnv &anticonn) override;
};


#endif //M7_FERMIONEXCITATIONGENERATOR_H
