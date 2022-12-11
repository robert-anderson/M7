//
// Created by Robert J. Anderson on 04/04/2022.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FrmExcitGen.h"

namespace exgen {
    struct UniformSingles : FrmExcitGen {
        UniformSingles(const FrmHam &h, PRNG &prng) :
                FrmExcitGen(h, prng, {opsig::c_sing}, "uniform") {}

        bool draw_frm(OpSig exsig, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) override;

        static bool draw_spin_conserve_fn(PRNG &prng, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn);

        static bool draw_spin_nonconserve_fn(PRNG &prng, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn);

        static prob_t prob_spin_conserve_fn(const field::FrmOnv &src, const conn::FrmOnv &conn);

        static prob_t prob_spin_nonconserve_fn(const field::FrmOnv &src, const conn::FrmOnv & /*conn*/);

        static prob_t prob_fn(const field::FrmOnv &src, const conn::FrmOnv &conn, bool spin_conserve);

    public:

        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;

        uint_t approx_nconn(OpSig /*exsig*/, sys::Particles particles) const override;
    };
}

#endif //M7_UNIFORMSINGLES_H
