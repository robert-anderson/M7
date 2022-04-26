//
// Created by anderson on 04/04/2022.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FrmExcitGen.h"

struct UniformSingles : FrmExcitGen {
    UniformSingles(const FrmHam &h, size_t nelec, PRNG &prng):
            FrmExcitGen(h, nelec, prng, {exsig_utils::ex_single}, "uniform"){}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

    static bool draw_spin_conserve_fn(PRNG &prng, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);

    static bool draw_spin_nonconserve_fn(PRNG &prng, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);
};


#endif //M7_UNIFORMSINGLES_H
