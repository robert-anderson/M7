//
// Created by rja on 04/08/2021.
//

#ifndef M7_FRMBOSHOLSTEIN_H
#define M7_FRMBOSHOLSTEIN_H

#include "ExcitGen.h"

struct FrmBosHolstein : public FrmBosExcitGen {
    const bool m_cre;
public:
    FrmBosHolstein(const Hamiltonian &h, PRNG &prng, bool cre) :
        FrmBosExcitGen(h, prng, {exsig_utils::ex_0001, exsig_utils::ex_0010}), m_cre(cre) {}

    bool draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override {

        if(!m_h.m_nboson_max) return false;

        const auto& occs = orbs.occ(src.m_frm).m_flat;
        DEBUG_ASSERT_EQ(src.m_bos.nelement(), src.m_frm.m_nsite,
                        "excit gen assumes one boson mode per fermion site");

        auto imode = occs[m_prng.draw_uint(occs.size())] % src.m_frm.m_nsite;
        // if m_cre, attempt to generate an ONV with an additional boson occupying the mode at imode
        // else, attempt to generate a "de-excited" ONV with one less boson occupying the mode at imode
        size_t curr_occ = src.m_bos[imode];
        DEBUG_ASSERT_LE(curr_occ, m_h.m_nboson_max, "current occupation of selected mode exceeds cutoff");

        // doubly-occupied sites are twice as likely to be drawn
        prob = src.m_frm.site_nocc(imode);
        prob/=occs.size();

        if (m_cre && curr_occ == m_h.m_nboson_max) return false;
        if (!m_cre && curr_occ == 0) return false;

        DEBUG_ASSERT_LE(size_t(curr_occ+m_cre), m_h.m_nboson_max, "generated boson occupation exceeds cutoff");

        conn.clear();
        if (m_cre) conn.m_bos.m_cre.add({imode, 1ul});
        else conn.m_bos.m_ann.add({imode, 1ul});
        return true;
    }

private:
    size_t approx_nconn() const override;
};


#endif //M7_FRMBOSHOLSTEIN_H
