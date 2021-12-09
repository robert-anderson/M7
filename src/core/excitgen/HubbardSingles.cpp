//
// Created by rja on 02/05/2021.
//

#include "HubbardSingles.h"

size_t HubbardSingles::approx_nconn() const {
    return m_nelec;
}

bool HubbardSingles::draw_frm(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs,
                            defs::prob_t &prob, conn::FrmOnv &conn) {
    auto rand = m_prng.draw_uint(m_nelec);
    const auto occ = orbs.occ(src).m_flat[rand];
    const auto isite = src.isite(occ);
    const auto ispin = src.ispin(occ);
    auto t_mat_row = h_cast()->m_t_mat_sparse[isite];
    const auto nvac = t_mat_row.first.size();
    rand = m_prng.draw_uint(nvac);
    auto vac = src.m_format.flatten({ispin, t_mat_row.first[rand]});
    if (src.get(vac)) return false;
    prob = 1.0 / double (m_nelec * nvac);
    conn.set(occ, vac);
    return true;
}

bool HubbardSingles::draw_frmbos(const size_t &exsig, const FrmBosOnv &src_onv, CachedOrbs &orbs,
                            defs::prob_t &prob, conn::FrmBosOnv &conn) {
    return draw(exsig, src_onv.m_frm, orbs, prob, conn.m_frm);
}

std::string HubbardSingles::description() const {
    return "hubbard hopping";
}
