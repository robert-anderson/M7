//
// Created by Robert J. Anderson on 05/04/2022.
//

#ifndef M7_BOSEXCITGEN_H
#define M7_BOSEXCITGEN_H

#include "M7_lib/hamiltonian/bos/BosHam.h"
#include "M7_lib/excitgen/ExcitGen.h"

#include <utility>

struct BosExcitGen : ExcitGen {

    const BosHam &m_h;

    BosExcitGen(const BosHam &h, PRNG &prng, defs::inds_t exsigs, std::string description);

    bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                     defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob,
                       defs::ham_t &helem, conn::FrmBosOnv &conn) override;

    bool draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::BosOnv &conn) override;

};


#endif //M7_BOSEXCITGEN_H