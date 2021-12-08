//
// Created by rja on 30/08/2021.
//

#include "LadderPureHolsteinZpm.h"

bool LadderPureHolsteinZpm::draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                 conn::FrmBosOnv &conn) {
    const auto& sites = orbs.not_singly_occ_sites(src.m_frm);
    if (sites.empty()) return false;
    size_t imode = sites[m_prng.draw_uint(sites.size())];

    bool cre = exsig_utils::decode_nbos_cre(exsig);

    size_t curr_occ = src.m_bos[imode];
    if (cre && curr_occ == m_h.m_nboson_max) return false;
    if (!cre && curr_occ == 0) return false;

    DEBUG_ASSERT_LE(size_t(curr_occ + cre), m_h.m_nboson_max, "generated boson occupation exceeds cutoff");

    conn.clear();
    if (cre) conn.m_bos.m_cre.add({imode, 1ul});
    else conn.m_bos.m_ann.add({imode, 1ul});
    prob = 1.0/sites.size();
    return true;
}

std::string LadderPureHolsteinZpm::description() const {
    return "ZPM-removed Holstein half-filling";
}
