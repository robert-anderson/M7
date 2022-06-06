//
// Created by Robert J. Anderson on 05/04/2022.
//

#include "HeisenbergExchange.h"

HeisenbergExchange::HeisenbergExchange(const FrmHam &h, PRNG &prng) :
        FrmLatticeExcitGen(h, prng, {exsig_utils::ex_double}, "lattice local exchange"){
    REQUIRE_TRUE(h.is<HeisenbergFrmHam>(), "given hamiltonian is not of HeisenbergFrmHam type");
}

bool HeisenbergExchange::draw_frm(const size_t &exsig, const field::FrmOnv &src,
                                  defs::prob_t &prob, conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig_utils::ex_double, "this excitation generator is only suitable for exsig 2200");
    const auto& lattice = m_h.m_basis.m_lattice;
    /*
     * the number of neighboring sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto& nconn_product = lattice->m_unique_nadj_product;
    const auto rand = m_prng.draw_uint(m_h.m_basis.m_nsite*nconn_product);
    const auto isite = m_h.m_basis.isite(rand / nconn_product);
    const auto ispin = src.get({1, isite});
    lattice->get_adj_row(isite, m_work_adj_row);
    const auto nvac = m_work_adj_row.size();
    const auto adj_elem = m_work_adj_row[rand % nvac];
    const auto jsite = adj_elem.m_isite;
    const auto jspin = src.get({1, jsite});
    if (jspin == ispin) return false; // no exchange
    DEBUG_ASSERT_NE(src.get({0, isite}), src.get({0, jsite}), "sites do not have opposite spins");
    // numerator of 2 since the exchange can be drawn in two ways
    prob = 2.0 / double (m_h.m_basis.m_nsite * nvac);
    //    <i  j  |  a  b>
    const auto i = m_h.m_basis.ispinorb(ispin, isite);
    const auto j = m_h.m_basis.ispinorb(jspin, jsite);
    const auto a = m_h.m_basis.ispinorb(ispin, jsite);
    const auto b = m_h.m_basis.ispinorb(jspin, isite);
    conn.m_ann.set_in_order(i, j);
    conn.m_cre.set_in_order(a, b);
    return true;
}

size_t HeisenbergExchange::approx_nconn(size_t exsig, sys::Particles particles) const {
    return 1ul<<m_h.m_basis.m_lattice->m_nadj_max;
}
