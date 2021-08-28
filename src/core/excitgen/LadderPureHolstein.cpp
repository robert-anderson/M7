//
// Created by rja on 04/08/2021.
//

#include "LadderPureHolstein.h"

bool LadderPureHolstein::draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                              conn::FrmBosOnv &conn) {

    if (!m_h.m_nboson_max) return false;

    bool cre = exsig_utils::decode_nbos_cre(exsig);

    const auto &occs = orbs.occ(src.m_frm).m_flat;
    DEBUG_ASSERT_EQ(src.m_bos.nelement(), src.m_frm.m_nsite,
                    "excit gen assumes one boson mode per fermion site");

    auto imode = occs[m_prng.draw_uint(occs.size())] % src.m_frm.m_nsite;
    // if m_cre, attempt to generate an ONV with an additional boson occupying the mode at imode
    // else, attempt to generate a "de-excited" ONV with one less boson occupying the mode at imode
    size_t curr_occ = src.m_bos[imode];
    DEBUG_ASSERT_LE(curr_occ, m_h.m_nboson_max, "current occupation of selected mode exceeds cutoff");

    // doubly-occupied sites are twice as likely to be drawn
    prob = src.m_frm.site_nocc(imode);
    prob /= occs.size();

    if (cre && curr_occ == m_h.m_nboson_max) return false;
    if (!cre && curr_occ == 0) return false;

    DEBUG_ASSERT_LE(size_t(curr_occ + cre), m_h.m_nboson_max, "generated boson occupation exceeds cutoff");

    conn.clear();
    if (cre) conn.m_bos.m_cre.add({imode, 1ul});
    else conn.m_bos.m_ann.add({imode, 1ul});
    return true;
}
