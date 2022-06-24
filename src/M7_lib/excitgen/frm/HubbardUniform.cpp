//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "HubbardUniform.h"

HubbardUniform::HubbardUniform(const FrmHam &h, PRNG &prng) :
        FrmLatticeExcitGen(h, prng, {exsig::ex_single}, "hubbard hopping") {
    REQUIRE_TRUE(h.is<HubbardFrmHam>(), "given hamiltonian is not of HubbardFrmHam type");
}

bool HubbardUniform::draw_frm(uint_t, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) {
    const auto& h = *m_h.as<HubbardFrmHam>();
    /*
     * the number of adjacent sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto nconn_product = h.m_basis.m_lattice->m_unique_nadj_product;
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto nelec = occs.size();
    /*
     * draw a random number with enough entropy to choose both the occupied spin orbital and a connected site regardless
     * of the actual connectivity
     */
    auto rand = m_prng.draw_uint(nelec * nconn_product);
    const auto occ = occs[rand / nconn_product];
    const auto isite = src.m_basis.isite(occ);
    const auto ispin = src.m_basis.ispin(occ);
    /*
     * fill the working vector with the adjacency of the occupied site
     */
    h.m_basis.m_lattice->get_adj_row(isite, m_work_adj_row);
    const auto nvac = m_work_adj_row.size();
    DEBUG_ASSERT_TRUE(nvac, "a lattice site should always have a non-zero number of adjacent sites");
    auto vac = src.m_format.flatten({ispin, m_work_adj_row[rand % nvac].m_isite});
    /*
     * return a null excitation if the drawn adjacent site is occupied in the drawn spin
     */
    if (src.get(vac)) return false;
    prob = 1.0 / double(nelec * nvac);
    conn.m_ann.set(occ);
    conn.m_cre.set(vac);
    return true;
}

prob_t HubbardUniform::prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const {
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto nelec = occs.size();
    const auto isite = src.m_basis.isite(conn.m_ann[0]);
    return 1.0/(src.m_basis.m_lattice->m_nadjs[isite]*nelec);
}

uint_t HubbardUniform::approx_nconn(uint_t, sys::Particles particles) const {
    return particles.m_frm;
}