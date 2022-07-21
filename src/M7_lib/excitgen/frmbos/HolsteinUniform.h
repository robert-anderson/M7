//
// Created by rja on 21/07/22.
//

#ifndef M7_EXGEN_HOLSTEIN_H
#define M7_EXGEN_HOLSTEIN_H

#include "FrmBosExcitGen.h"

namespace exgen {
    /**
     * Excitation generator for the Holstein fermion-boson Hamiltonian. Even though the exsigs drawn are 0010 and 0001,
     * BosOnvs are never successfully drawn since the matrix element is only non-zero with the mode is coordinated to an
     * occupied fermionic site
     */

    struct HolsteinUniform : FrmBosExcitGen {

        HolsteinUniform(const FrmBosHam& h, PRNG& prng, uint_t exsig);
        uint_t approx_nconn(uint_t /*exsig*/, sys::Particles particles) const override;

    };

    struct HolsteinUniform0010 : HolsteinUniform {

        HolsteinUniform0010(const FrmBosHam& h, PRNG& prng): HolsteinUniform(h, prng, exsig::ex_0010){}

        bool draw_frmbos(uint_t exsig, const field::FrmBosOnv& src, prob_t& prob, conn::FrmBosOnv& conn) override;

        prob_t prob_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn) const override;
    };

    struct HolsteinUniform0001 : HolsteinUniform {

        HolsteinUniform0001(const FrmBosHam& h, PRNG& prng): HolsteinUniform(h, prng, exsig::ex_0001){}

        bool draw_frmbos(uint_t exsig, const field::FrmBosOnv& src, prob_t& prob, conn::FrmBosOnv& conn) override;

        prob_t prob_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& /*conn*/) const override;
    };
}


#endif //M7_EXGEN_HOLSTEIN_H
