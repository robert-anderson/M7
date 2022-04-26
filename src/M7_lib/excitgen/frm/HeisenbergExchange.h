//
// Created by rja on 05/04/2022.
//

#ifndef M7_HEISENBERGEXCHANGE_H
#define M7_HEISENBERGEXCHANGE_H

#include "FrmExcitGen.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"

struct HeisenbergExchange : FrmExcitGen {

    HeisenbergExchange(const FrmHam &h, size_t nelec, PRNG &prng);

    virtual ~HeisenbergExchange() {}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

};


#endif //M7_HEISENBERGEXCHANGE_H
