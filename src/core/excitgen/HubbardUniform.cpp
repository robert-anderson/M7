//
// Created by rja on 02/05/2021.
//

#include "HubbardUniform.h"

size_t HubbardUniform::approx_nconn() const {
    return m_nelec;
}

bool HubbardUniform::draw_frm(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs,
                            defs::prob_t &prob, conn::FrmOnv &conn) {
    const auto h = h_cast();
    /*
     * the number of neighboring sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto& nconn_product = h->m_lattice.m_unique_nconn_product;
    auto rand = m_prng.draw_uint(m_nelec*nconn_product);
    const auto occ = orbs.occ(src).m_flat[rand/nconn_product];
    const auto isite = src.isite(occ);
    const auto ispin = src.ispin(occ);
    auto t_mat_row = h->m_lattice.m_sparse[isite];
    const auto nvac = t_mat_row.first.size();
    auto vac = src.m_format.flatten({ispin, t_mat_row.first[rand%nvac]});
    if (src.get(vac)) return false;
    prob = 1.0 / double (m_nelec * nvac);
    conn.set(occ, vac);
    return true;
}

bool HubbardUniform::draw_frmbos(const size_t &exsig, const FrmBosOnv &src_onv, CachedOrbs &orbs,
                            defs::prob_t &prob, conn::FrmBosOnv &conn) {
    return draw(exsig, src_onv.m_frm, orbs, prob, conn.m_frm);
}

std::string HubbardUniform::description() const {
    return "hubbard hopping";
}
