//
// Created by rja on 05/04/2022.
//

#ifndef M7_HEISENBERGEXCHANGE2_H
#define M7_HEISENBERGEXCHANGE2_H

#include "FrmExcitGen2.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"

struct HeisenbergExchange2 : FrmExcitGen2 {

    const HeisenbergFrmHam *h_cast() const;

    HeisenbergExchange2(const FrmHam &h, PRNG &prng);

    virtual ~HeisenbergExchange2() {}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

};


#endif //M7_HEISENBERGEXCHANGE2_H
