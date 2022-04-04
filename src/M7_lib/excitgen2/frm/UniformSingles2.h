//
// Created by anderson on 04/04/2022.
//

#ifndef M7_UNIFORMSINGLES2_H
#define M7_UNIFORMSINGLES2_H

#include "FrmExcitGen2.h"

struct UniformSingles2 : FrmExcitGen2 {
    UniformSingles2(const FrmHam &h, PRNG &prng):
        FrmExcitGen2(h, prng, {exsig_utils::ex_single}, "uniform"){}

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                    conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;
};


#endif //M7_UNIFORMSINGLES2_H
