//
// Created by Robert J. Anderson on 05/04/2022.
//

#ifndef M7_EXGEN_HEISENBERGEXCHANGE_H
#define M7_EXGEN_HEISENBERGEXCHANGE_H

#include "FrmExcitGen.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"

namespace exgen {
    struct HeisenbergExchange : FrmLatticeExcitGen {
        HeisenbergExchange(const FrmHam &h, PRNG &prng);

        virtual ~HeisenbergExchange() {}

        bool draw_frm(OpSig exsig, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) override;

        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;

        uint_t approx_nconn(OpSig /*exsig*/, sys::Particles /*particles*/) const override;

    };
}

#endif //M7_EXGEN_HEISENBERGEXCHANGE_H
