//
// Created by Robert J. Anderson on 04/04/2022.
//

#ifndef M7_HUBBARDUNIFORM_H
#define M7_HUBBARDUNIFORM_H

#include "FrmExcitGen.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"

/**
 * excitation generator for the N-dimensional Hubbard model which does not constrain the creation operator drawn based
 * on the occupation
 */
struct HubbardUniform : FrmLatticeExcitGen {

    HubbardUniform(const FrmHam& h, PRNG& prng);

    virtual ~HubbardUniform(){}

    bool draw_frm(uint_t /*exsig*/, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) override;

    prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;

    uint_t approx_nconn(uint_t exsig, sys::Particles particles) const override;

};


#endif //M7_HUBBARDUNIFORM_H
