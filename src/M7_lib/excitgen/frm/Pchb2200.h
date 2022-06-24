//
// Created by Robert J. Anderson on 04/04/2022.
//

//TODO: rename guard define
#ifndef M7_PCHBDOUBLES_H
#define M7_PCHBDOUBLES_H

#include "FrmExcitGen.h"
#include "M7_lib/sample/Aliaser.h"
#include "M7_lib/parallel/MPIAssert.h"

struct Pchb2200 : FrmExcitGen {
private:
    const uint_t m_nspinorb_pair;
    Aliaser m_pick_ab_given_ij;

public:
    Pchb2200(const FrmHam& h, PRNG& prng);

    bool draw_h_frm(uint_t exsig, const field::FrmOnv& src, prob_t& prob,
                    ham_t& helem, conn::FrmOnv& conn) override;

    bool draw_frm(uint_t exsig, const field::FrmOnv& src, prob_t& prob, conn::FrmOnv& conn) override;

    prob_t prob_h_frm(const field::FrmOnv& src, const conn::FrmOnv& conn, ham_t helem) const override;

    prob_t prob_frm(const field::FrmOnv& src, const conn::FrmOnv& conn) const override;

    uint_t approx_nconn(uint_t /*exsig*/, sys::Particles particles) const override;

};


#endif //M7_PCHBDOUBLES_H
