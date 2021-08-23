//
// Created by rja on 04/08/2021.
//

#ifndef M7_UNIFORMHOLSTEIN_H
#define M7_UNIFORMHOLSTEIN_H

#include "ExcitGen.h"

class UniformHolstein : public FrmBosExcitGen {
    const bool m_cre;
public:
    UniformHolstein(const Hamiltonian &h, PRNG &prng, bool cre) : FrmBosExcitGen(h, prng), m_cre(cre){}

    bool draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
              defs::ham_t &helem, conn::FrmBosOnv &conn) override;

private:
    size_t approx_nconn() const override;
};


#endif //M7_UNIFORMHOLSTEIN_H
