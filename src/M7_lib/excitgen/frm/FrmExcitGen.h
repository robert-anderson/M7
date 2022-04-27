//
// Created by Robert J. Anderson on 03/04/2022.
//

#ifndef M7_FRMEXCITGEN_H
#define M7_FRMEXCITGEN_H

#include "M7_lib/excitgen/ExcitGen.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct FrmExcitGen : ExcitGen {

    const FrmHam &m_h;
    const sys::frm::Sector m_sector;

    FrmExcitGen(const FrmHam &h, sys::frm::Electrons elecs, PRNG &prng, defs::inds exsigs, std::string description);

    bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                     defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob,
                       defs::ham_t &helem, conn::FrmBosOnv &conn) override;

    bool draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::BosOnv &conn) override;

};


#endif //M7_FRMEXCITGEN_H
