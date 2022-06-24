//
// Created by Robert J. Anderson on 05/04/2022.
//

#ifndef M7_HEISENBERGEXCHANGE_H
#define M7_HEISENBERGEXCHANGE_H

#include "FrmExcitGen.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"

struct HeisenbergExchange : FrmLatticeExcitGen {
    HeisenbergExchange(const FrmHam& h, PRNG& prng);

    virtual ~HeisenbergExchange() {}

    bool draw_frm(size_t exsig, const field::FrmOnv& src, defs::prob_t& prob, conn::FrmOnv& conn) override;

    defs::prob_t prob_frm(const field::FrmOnv& src, const conn::FrmOnv& conn) const override;

    size_t approx_nconn(size_t /*exsig*/, sys::Particles /*particles*/) const override;

};


#endif //M7_HEISENBERGEXCHANGE_H
