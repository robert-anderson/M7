//
// Created by rja on 04/08/2021.
//

#ifndef M7_UNIFORMFRMBOS_H
#define M7_UNIFORMFRMBOS_H

#include "ExcitGen.h"

class UniformFrmBos : public FrmBosExcitGen {
public:
    UniformFrmBos(const Hamiltonian &h, PRNG &prng) : FrmBosExcitGen(h, prng){}

    bool draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
              defs::ham_t &helem, conn::FrmBosOnv &conn) override;

private:
    std::string description() const override;

    size_t approx_nconn() const override;
};


#endif //M7_UNIFORMFRMBOS_H
