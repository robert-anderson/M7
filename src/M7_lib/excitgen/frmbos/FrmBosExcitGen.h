//
// Created by Robert J. Anderson on 05/04/2022.
//

#ifndef M7_FRMBOSEXCITGEN_H
#define M7_FRMBOSEXCITGEN_H

#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"
#include "M7_lib/excitgen/ExcitGen.h"

#include <utility>

struct FrmBosExcitGen : ExcitGen {

    const FrmBosHam& m_h;

    FrmBosExcitGen(const FrmBosHam& h, PRNG& prng, defs::uintv_t exsigs, std::string description);

    bool draw_h_frm(uint_t /*exsig*/, const field::FrmOnv& /*src*/, defs::prob_t& prob, defs::ham_t& helem,
                    conn::FrmOnv& /*conn*/) override {
        prob = 0.0;
        helem = 0.0;
        return false;
    }

    bool draw_h_frmbos(uint_t exsig, const field::FrmBosOnv& src, defs::prob_t& prob, defs::ham_t& helem,
                       conn::FrmBosOnv& conn) override;

    bool draw_h_bos(uint_t /*exsig*/, const field::BosOnv& /*src*/, defs::prob_t& prob, defs::ham_t& helem,
                    conn::BosOnv& /*conn*/) override {
        prob = 0.0;
        helem = 0.0;
        return false;
    }
};

#endif //M7_FRMBOSEXCITGEN_H