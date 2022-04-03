//
// Created by rja on 04/04/2022.
//

#ifndef M7_PCHBDOUBLES_H
#define M7_PCHBDOUBLES_H

#include "FrmExcitGen2.h"
#include "M7_lib/sample/Aliaser.h"
#include "M7_lib/parallel/MPIAssert.h"

struct PchbDoubles : FrmExcitGen2 {
private:
    Aliaser m_pick_ab_given_ij;

public:
    PchbDoubles(const FrmHam &h, PRNG &prng);

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

};


#endif //M7_PCHBDOUBLES_H
