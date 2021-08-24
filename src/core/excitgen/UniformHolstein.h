//
// Created by rja on 04/08/2021.
//

#ifndef M7_UNIFORMHOLSTEIN_H
#define M7_UNIFORMHOLSTEIN_H

#include "ExcitGen.h"

class UniformHolstein : public FrmBosExcitGen {
    const bool m_cre;
public:
    UniformHolstein(const Hamiltonian &h, PRNG &prng, bool cre) : FrmBosExcitGen(h, prng), m_cre(cre) {}

    bool draw(const size_t &exsig, const FrmBosOnv &onv, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override {

        if(!m_h.m_nboson_max) return false;

        const auto& occs = orbs.occ(onv.m_frm).m_flat;
        DEBUG_ASSERT_EQ(onv.m_bos.nelement(), onv.m_frm.m_nsite,
                        "excit gen assumes one boson mode per fermion site");

        auto imode = occs[m_prng.draw_uint(occs.size())] % onv.m_frm.m_nsite;
        // if m_cre, attempt to generate an ONV with an additional boson occupying the mode at imode
        // else, attempt to generate a "de-excited" ONV with one less boson occupying the mode at imode
        size_t curr_occ = onv.m_bos[imode];
        DEBUG_ASSERT_LE(curr_occ, m_h.m_nboson_max, "current occupation of selected mode exceeds cutoff");

        // doubly-occupied sites are twice as likely to be drawn
        prob = onv.m_frm.site_nocc(imode);
        prob/=occs.size();

        if (m_cre && curr_occ == m_h.m_nboson_max) return false;
        if (!m_cre && curr_occ == 0) return false;

        DEBUG_ASSERT_LE(size_t(curr_occ+m_cre), m_h.m_nboson_max, "generated boson occupation exceeds cutoff");

        conn.clear();
        if (m_cre) conn.m_bos.m_cre.add({imode, 1ul});
        else conn.m_bos.m_ann.add({imode, 1ul});

        helem = m_h.get_element(onv, conn);
        return true;
    }

private:
    size_t approx_nconn() const override;
};


#endif //M7_UNIFORMHOLSTEIN_H
