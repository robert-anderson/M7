//
// Created by anderson on 2/6/22.
//

#include "HeisenbergUniform.h"

bool HeisenbergUniform::draw_frm(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                 conn::FrmOnv &conn) {
    const auto h = h_cast();
    /*
     * the number of neighboring sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto& nconn_product = h->m_lattice.m_unique_nconn_product;
    auto rand = m_prng.draw_uint(h->m_nsite*nconn_product);
    const auto isite = src.isite(rand/nconn_product);
    const auto ispin = src.ispin(isite);
    auto t_mat_row = h->m_lattice.m_sparse[isite];
    const auto nvac = t_mat_row.first.size();
    auto jsite = t_mat_row.first[rand%nvac];
    if (src.ispin(jsite)==ispin) return false; // no exchange
    prob = 1.0 / double (h->m_nsite * nvac);
    conn.set(isite, jsite, jsite, isite);
    return true;
}

std::string HeisenbergUniform::description() const {
    return "heisenberg exchange";
}
