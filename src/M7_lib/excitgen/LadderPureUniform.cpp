//
// Created by rja on 28/08/2021.
//

#include "LadderPureUniform.h"

bool LadderPureUniform::draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                    conn::FrmBosOnv &conn) {
    if (!m_h.m_nboson_max) return false;

    bool cre = exsig_utils::decode_nbos_cre(exsig);
    auto imode = m_prng.draw_uint(src.m_bos.m_nelement);

    size_t curr_occ = src.m_bos[imode];
    if (cre && curr_occ == m_h.m_nboson_max) return false;
    if (!cre && curr_occ == 0) return false;

    DEBUG_ASSERT_LE(size_t(curr_occ + cre), m_h.m_nboson_max, "generated boson occupation exceeds cutoff");

    conn.clear();
    if (cre) conn.m_bos.m_cre.set(imode);
    else conn.m_bos.m_ann.set(imode);
    prob = 1.0/src.m_bos.m_nelement;
    return true;
}


size_t LadderPureUniform::approx_nconn() const {
    // assume there's one excitation or de-excitation available per electron
    return m_h.m_frm->m_nelec;
}

std::string LadderPureUniform::description() const {
    return "uniform";
}
