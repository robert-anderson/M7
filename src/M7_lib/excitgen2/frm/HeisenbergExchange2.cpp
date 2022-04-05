//
// Created by rja on 05/04/2022.
//

#include "HeisenbergExchange2.h"

const HeisenbergFrmHam *HeisenbergExchange2::h_cast() const {
    return dynamic_cast<const HeisenbergFrmHam *>(&m_h);
}

HeisenbergExchange2::HeisenbergExchange2(const FrmHam &h, PRNG &prng) :
        FrmExcitGen2(h, prng, {exsig_utils::ex_double}, "lattice local exchange"){}

bool HeisenbergExchange2::draw_frm(const size_t &exsig, const field::FrmOnv &src,
                                   defs::prob_t &prob, conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig_utils::ex_double, "this excitation generator is only suitable for exsig 2200");
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
    const auto& row = h->m_lattice.m_sparse[isite];
    const auto nvac = row.first.size();
    const auto jsite = row.first[rand % nvac];
    const auto jspin = src.get({1, jsite});
    if (jspin == ispin) return false; // no exchange
    DEBUG_ASSERT_NE(src.get({0, isite}), src.get({0, jsite}), "sites do not have opposite spins");
    // numerator of 2 since the exchange can be drawn in two ways
    prob = 2.0 / double (h->m_nsite * nvac);
    //    <i  j  |  a  b>
    const auto i = src.ibit(ispin, isite);
    const auto j = src.ibit(jspin, jsite);
    const auto a = src.ibit(ispin, jsite);
    const auto b = src.ibit(jspin, isite);
    if (jspin) conn.set(i, j, a, b);
    else conn.set(j, i, b, a);
    return true;
}

size_t HeisenbergExchange2::approx_nconn() const {
    return 1ul<<h_cast()->m_lattice.m_spec.m_format.m_nind;
}
