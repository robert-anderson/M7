//
// Created by rja on 04/04/2022.
//

#include "HubbardUniform2.h"

HubbardUniform2::HubbardUniform2(const FrmHam &h, PRNG &prng) :
    FrmExcitGen2(h, prng, {exsig_utils::ex_single}, "hubbard hopping"){
    REQUIRE_TRUE(h_cast(), "given hamiltonian is not of HubbardFrmHam type");
}

bool HubbardUniform2::draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
    const auto h = h_cast();
    /*
     * the number of neighboring sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto& nconn_product = h->m_lattice.m_unique_nconn_product;
    auto rand = m_prng.draw_uint(h->m_nelec*nconn_product);
    const auto occ = src.m_decoded.m_simple_occs.get()[rand/nconn_product];
    const auto isite = src.isite(occ);
    const auto ispin = src.ispin(occ);
    auto t_mat_row = h->m_lattice.m_sparse[isite];
    const auto nvac = t_mat_row.first.size();
    auto vac = src.m_format.flatten({ispin, t_mat_row.first[rand%nvac]});
    if (src.get(vac)) return false;
    prob = 1.0 / double (h->m_nelec * nvac);
    conn.set(occ, vac);
    return true;
}

size_t HubbardUniform2::approx_nconn() const {
    return m_h.m_nelec;
}