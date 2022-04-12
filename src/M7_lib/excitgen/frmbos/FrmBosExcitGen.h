//
// Created by rja on 05/04/2022.
//

#ifndef M7_FRMBOSEXCITGEN_H
#define M7_FRMBOSEXCITGEN_H

#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"
#include "M7_lib/excitgen/ExcitGen.h"

#include <utility>

struct FrmBosExcitGen : ExcitGen {

    const FrmBosHam& m_h;

    FrmBosExcitGen(const FrmBosHam& h, PRNG &prng, defs::inds exsigs, std::string description);

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                    conn::FrmOnv &conn) override;

    bool draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                       conn::FrmBosOnv &conn) override;

    bool draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                    conn::BosOnv &conn) override;
};


#endif //M7_FRMBOSEXCITGEN_H
