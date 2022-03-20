//
// Created by rja on 04/08/2021.
//

#include "LadderPureHolstein.h"

std::string LadderPureHolstein::description() const {
    return "holstein";
}

bool LadderHolsteinCre::draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                     conn::FrmBosOnv &conn) {
    if (!m_h.m_nboson_max) return false;
    const auto &occs = orbs.occ(src.m_frm).m_flat.inds();

    auto imode = src.m_frm.isite(occs[m_prng.draw_uint(occs.size())]);
    // attempt to generate an ONV with an additional boson occupying the mode at imode
    size_t curr_occ = src.m_bos[imode];
    DEBUG_ASSERT_LE(curr_occ, m_h.m_nboson_max, "current occupation of selected mode exceeds cutoff");

    // doubly-occupied sites are twice as likely to be drawn
    prob = src.m_frm.site_nocc(imode);
    prob /= occs.size();

    if (curr_occ == m_h.m_nboson_max) return false;

    conn.clear();
    conn.m_bos.m_cre.add({imode, 1ul});
    return true;
}

bool LadderHolsteinAnn::draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                    conn::FrmBosOnv &conn) {
    if (!m_h.m_nboson_max) return false;

    const auto &occs = orbs.occ_sites_nonzero_bosons(src);
    if (occs.empty()) return false;

    auto imode = occs[m_prng.draw_uint(occs.size())];
    // attempt to generate a "de-excited" ONV with one less boson occupying the mode at imode
    DEBUG_ASSERT_LE(src.m_bos[imode], m_h.m_nboson_max, "current occupation of selected mode exceeds cutoff");
    DEBUG_ASSERT_TRUE(src.m_bos[imode], "selected boson mode should have non-zero occupation");

    prob = 1.0/occs.size();

    conn.clear();
    conn.m_bos.m_ann.add({imode, 1ul});
    return true;
}
