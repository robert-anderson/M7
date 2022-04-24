//
// Created by rja on 04/04/2022.
//

#ifndef M7_HUBBARDUNIFORM_H
#define M7_HUBBARDUNIFORM_H

#include "FrmExcitGen.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"

struct HubbardUniform : FrmExcitGen {

    HubbardUniform(const FrmHam& h, PRNG& prng);

    virtual ~HubbardUniform(){}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

};


#endif //M7_HUBBARDUNIFORM_H
