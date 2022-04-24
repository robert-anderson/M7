//
// Created by rja on 05/04/2022.
//

#include "HeisenbergExchange.h"

HeisenbergExchange::HeisenbergExchange(const FrmHam &h, PRNG &prng) :
        FrmExcitGen(h, prng, {exsig_utils::ex_double}, "lattice local exchange"){
    REQUIRE_TRUE(h.is<HeisenbergFrmHam>(), "given hamiltonian is not of HeisenbergFrmHam type");
}

bool HeisenbergExchange::draw_frm(const size_t &exsig, const field::FrmOnv &src,
                                  defs::prob_t &prob, conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig_utils::ex_double, "this excitation generator is only suitable for exsig 2200");
    const auto& lattice = m_h.as<HeisenbergFrmHam>()->m_lattice;
    /*
     * the number of neighboring sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto& nconn_product = lattice.m_unique_nconn_product;
    const auto rand = m_prng.draw_uint(m_h.m_hs.m_sites*nconn_product);
    const auto isite = m_h.m_hs.m_sites.isite(rand / nconn_product);
    const auto ispin = src.get({1, isite});
    const auto& row = lattice.m_sparse[isite];
    const auto nvac = row.first.size();
    const auto jsite = row.first[rand % nvac];
    const auto jspin = src.get({1, jsite});
    if (jspin == ispin) return false; // no exchange
    DEBUG_ASSERT_NE(src.get({0, isite}), src.get({0, jsite}), "sites do not have opposite spins");
    // numerator of 2 since the exchange can be drawn in two ways
    prob = 2.0 / double (m_h.m_hs.m_sites * nvac);
    //    <i  j  |  a  b>
    const auto i = m_h.m_hs.m_sites.ispinorb(ispin, isite);
    const auto j = m_h.m_hs.m_sites.ispinorb(jspin, jsite);
    const auto a = m_h.m_hs.m_sites.ispinorb(ispin, jsite);
    const auto b = m_h.m_hs.m_sites.ispinorb(jspin, isite);
    conn.m_ann.set_in_order(i, j);
    conn.m_cre.set_in_order(a, b);
    return true;
}

size_t HeisenbergExchange::approx_nconn() const {
    return 1ul<<m_h.as<HeisenbergFrmHam>()->m_lattice.m_spec.m_format.m_nind;
}
