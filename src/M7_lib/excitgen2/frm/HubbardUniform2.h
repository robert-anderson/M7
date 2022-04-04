//
// Created by rja on 04/04/2022.
//

#ifndef M7_HUBBARDUNIFORM2_H
#define M7_HUBBARDUNIFORM2_H

#include "FrmExcitGen2.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"

struct HubbardUniform2 : FrmExcitGen2 {

    const HubbardFrmHam* h_cast() const {
        return dynamic_cast<const HubbardFrmHam*>(&m_h);
    }

    HubbardUniform2(const FrmHam& h, PRNG& prng);

    virtual ~HubbardUniform2(){}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

};


#endif //M7_HUBBARDUNIFORM2_H
