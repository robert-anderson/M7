//
// Created by Robert J. Anderson on 04/04/2022.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FrmExcitGen.h"

struct UniformSingles : FrmExcitGen {
    UniformSingles(const FrmHam &h, PRNG &prng):
            FrmExcitGen(h, prng, {exsig_utils::ex_single}, "uniform"){}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    static bool draw_spin_conserve_fn(PRNG &prng, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);

    static bool draw_spin_nonconserve_fn(PRNG &prng, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);

    size_t approx_nconn(size_t exsig, sys::Particles particles) const override;
};


#endif //M7_UNIFORMSINGLES_H
