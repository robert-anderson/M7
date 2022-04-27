//
// Created by Robert J. Anderson on 2/6/22.
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
    const auto rand = m_prng.draw_uint(h->m_nsite*nconn_product);
    const auto isite = src.isite(rand / nconn_product);
    const auto ispin = src.get({1, isite});
    const auto & row = h->m_lattice.m_sparse[isite];
    const auto nvac = row.first.size();
    const auto jsite = row.first[rand % nvac];
    const auto jspin = src.get({1, jsite});
    if (jspin == ispin) return false; // no exchange
    DEBUG_ASSERT_NE(src.get({0, isite}), src.get({0, jsite}), "sites do not have opposite spins");
    prob = 1.0 / double (h->m_nsite * nvac);
    //    <i  j  |  a  b>
    const auto i = src.ibit(ispin, isite);
    const auto j = src.ibit(jspin, jsite);
    const auto a = src.ibit(ispin, jsite);
    const auto b = src.ibit(jspin, isite);
    conn.m_ann.set_in_order(i, j);
    conn.m_cre.set_in_order(a, b);
    return true;
}

std::string HeisenbergUniform::description() const {
    return "heisenberg exchange";
}

size_t HeisenbergUniform::approx_nconn() const {
    return 1ul<<h_cast()->m_lattice.m_spec.m_format.m_nind;
}
