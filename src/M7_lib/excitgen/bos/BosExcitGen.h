//
// Created by Robert J. Anderson on 05/04/2022.
//

#ifndef M7_BOSEXCITGEN_H
#define M7_BOSEXCITGEN_H

#include "M7_lib/hamiltonian/bos/BosHam.h"
#include "M7_lib/excitgen/ExcitGen.h"

#include <utility>

struct BosExcitGen : ExcitGen {

    const BosHam& m_h;

    BosExcitGen(const BosHam& h, PRNG& prng, uintv_t exsigs, str_t description);

    bool draw_frmbos(uint_t exsig, const field::FrmBosOnv& src,
                     prob_t& prob, conn::FrmBosOnv& conn) override;

    bool draw_h_frm(uint_t /*exsig*/, const field::FrmOnv& /*src*/, prob_t& prob,
                    ham_t& helem, conn::FrmOnv& /*conn*/) override {
        prob = 0.0;
        helem = 0.0;
        return false;
    }

    bool draw_h_frmbos(uint_t exsig, const field::FrmBosOnv& src, prob_t& prob,
                       ham_t& helem, conn::FrmBosOnv& conn) override;

    bool draw_h_bos(uint_t exsig, const field::BosOnv& src, prob_t& prob,
                    ham_t& helem, conn::BosOnv& conn) override;

};


#endif //M7_BOSEXCITGEN_H
