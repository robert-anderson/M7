//
// Created by rja on 03/04/2022.
//

#ifndef M7_FRMEXCITGEN2_H
#define M7_FRMEXCITGEN2_H

#include "M7_lib/excitgen2/ExcitGen2.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct FrmExcitGen2 : ExcitGen2 {

    const FrmHam &m_h;
    /**
     * number of pairs of electrons in the system (FrmHam conserves fermion number)
     */
    const size_t m_nelec_pair;
    /**
     * number of pairs of fermionic degrees of freedom in the system
     */
    const size_t m_nspinorb_pair;

    FrmExcitGen2(const FrmHam &h, PRNG &prng, defs::inds exsigs, std::string description);

    bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                     defs::prob_t &prob, conn::FrmBosOnv &conn) override;

    bool draw_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) override;

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob,
                       defs::ham_t &helem, conn::FrmBosOnv &conn) override;

    bool draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob,
                    defs::ham_t &helem, conn::BosOnv &conn) override;

};


#endif //M7_FRMEXCITGEN2_H
