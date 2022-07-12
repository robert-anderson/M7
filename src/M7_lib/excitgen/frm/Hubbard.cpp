//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "Hubbard.h"

exgen::HubbardBase::HubbardBase(const FrmHam &h, PRNG &prng, str_t description) :
        FrmLatticeExcitGen(h, prng, {exsig::ex_single}, description) {
    REQUIRE_TRUE(h.is<HubbardFrmHam>(), "given hamiltonian is not of HubbardFrmHam type");
}

bool exgen::HubbardBase::draw_frm(uint_t exsig, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) {
    DEBUG_ONLY(exsig);
    DEBUG_ASSERT_EQ(exsig, exsig::ex_single, "this excitation generator is only suitable for exsig 1100");
    const auto& h = *m_h.as<HubbardFrmHam>();
    /*
     * the number of adjacent sites accessible is not decided till the occupied index has been chosen. If the integer
     * picked is an integral multiple of all possible numbers of accessible sites, then in any case the modular
     * remainder will provide an unbiased index - saving a PRNG call
     */
    const auto nconn_lcm = h.m_basis.m_lattice->m_lcm_le_nadj_max;
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto nelec = occs.size();
    /*
     * draw a random number with enough entropy to choose both the occupied spin orbital and a connected site regardless
     * of the actual connectivity
     */
    auto rand = m_prng.draw_uint(nelec * nconn_lcm);
    const auto occ = get_occ(rand, src, prob);
    DEBUG_ASSERT_GE(prob, 0.0, "prob is non-positive");
    DEBUG_ASSERT_LE(prob, 1.0, "prob is more than 1");
    const auto isite = src.m_basis.isite(occ);
    const auto ispin = src.m_basis.ispin(occ);
    /*
     * fill the working vector with the adjacency of the occupied site
     */
    h.m_basis.m_lattice->get_adj_row(isite, m_work_adj_row);
    /*
     * when selecting the vacant site, skip the sites which are occupied in the same spin channel as the chosen electron
     */
    set_valid_adj_vacant(src, ispin);
    const auto nvac = m_valid_in_adj_row.size();
    /*
     * if there are no vacants from which to choose, this occupied pick is a dead end
     */
    if (!nvac) return false;
    auto ivalid = rand % nvac;
    auto vac = src.m_format.flatten({ispin, m_valid_in_adj_row[ivalid]->m_isite});
    DEBUG_ASSERT_FALSE(src.get(vac), "should not have picked an occupied spin orb for the vacant");
    prob /= double(nvac);
    conn.m_ann.set(occ);
    conn.m_cre.set(vac);
    return true;
}

uint_t exgen::HubbardBase::approx_nconn(uint_t, sys::Particles particles) const {
    return particles.m_frm;
}