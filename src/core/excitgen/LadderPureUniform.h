//
// Created by rja on 28/08/2021.
//

#ifndef M7_LADDERPUREUNIFORM_H
#define M7_LADDERPUREUNIFORM_H

#include "ExcitGen.h"

struct LadderPureUniform : public LadderExcitGen {
    LadderPureUniform(const Hamiltonian& ham, PRNG& prng):
        LadderExcitGen(ham, prng, {exsig_utils::ex_0010, exsig_utils::ex_0001}){}

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    std::string description() const override;

private:
    size_t approx_nconn() const override;
};


#endif //M7_LADDERPUREUNIFORM_H
