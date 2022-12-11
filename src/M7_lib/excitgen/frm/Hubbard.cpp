//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "Hubbard.h"

exgen::HubbardBase::HubbardBase(const FrmHam &h, PRNG &prng, str_t description) :
        FrmLatticeExcitGen(h, prng, {opsig::c_sing}, description) {
    REQUIRE_TRUE(h.is<HubbardFrmHam>(), "given hamiltonian is not of HubbardFrmHam type");
}

uint_t exgen::HubbardBase::get_occ_uniform(uint32_t rand, const field::FrmOnv& src, prob_t& prob) const {
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto nconn_lcm = m_h.m_basis.m_lattice->m_lcm_le_nadj_max;
    const auto nelec = occs.size();
    prob = 1/double(nelec);
    return occs[rand / nconn_lcm];
}

prob_t exgen::HubbardBase::prob_uniform(const field::FrmOnv& src, const conn::FrmOnv& conn) const {
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto nelec = occs.size();
    const auto isite = src.m_basis.isite(conn.m_ann[0]);
    const auto ispin = src.m_basis.ispin(conn.m_ann[0]);

    auto& valid_adj = this->valid_adj(isite, src, ispin);
    return 1.0/double(valid_adj.size()*nelec);
}


bool exgen::HubbardBase::draw_frm(OpSig exsig, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) {
    DEBUG_ONLY(exsig);
    DEBUG_ASSERT_EQ(exsig, opsig::c_sing, "this excitation generator is only suitable for exsig 1100");
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

    auto& valid_adj = this->valid_adj(isite, src, ispin);
    const auto nvac = valid_adj.size();
    /*
     * if there are no vacants from which to choose, this occupied pick is a dead end
     */
    if (!nvac) return false;
    auto ivalid = rand % nvac;
    auto vac = src.m_format.flatten({ispin, valid_adj[ivalid]->m_i});
    DEBUG_ASSERT_FALSE(src.get(vac), "should not have picked an occupied spin orb for the vacant");
    prob /= double(nvac);
    conn.m_ann.set(occ);
    conn.m_cre.set(vac);
    return true;
}

uint_t exgen::HubbardBase::approx_nconn(OpSig, sys::Particles particles) const {
    return particles.m_frm;
}

exgen::HubbardPreferDoubleOcc::HubbardPreferDoubleOcc(const FrmHam& h, PRNG& prng, prob_t doub_occ_u_fac) :
        HubbardBase(h, prng, "doubly occupied site-preferring hubbard hopping"),
        m_prob_doub_occ(1.0/(1.0+1.0/(doub_occ_u_fac*m_h.as<HubbardFrmHam>()->m_u))){}

prob_t exgen::HubbardPreferDoubleOcc::combined_occ_prob(uint_t nelec, uint_t nelec_doub_occ) const {
    // prob = prob_uni + prob_pref
    return (1.0 - m_prob_doub_occ)/prob_t(nelec) + m_prob_doub_occ/prob_t(nelec_doub_occ);
}

uint_t exgen::HubbardPreferDoubleOcc::get_occ(uint32_t rand, const field::FrmOnv& src, prob_t& prob) const {
    /*
     * total number of electrons
     */
    const auto nelec = src.m_decoded.m_simple_occs.get().size();
    const auto& doub_occs = src.m_decoded.m_doubly_occ_sites.get();
    /*
     * if there are no double occupations, then do uniform selection without alteration of prob
     */
    if (doub_occs.empty()) return get_occ_uniform(rand, src, prob);
    /*
     * number of electrons in doubly occupied sites
     */
    const auto nelec_doub_occ = 2*doub_occs.size();
    /*
     * it should still be possible to draw any electron
     */
    if (m_prng.draw_float() > m_prob_doub_occ) {
        /*
         * uniform branch: draw using the uniform strategy
         */
        const auto occ = get_occ_uniform(rand, src, prob);
        const auto isite = src.m_basis.isite(occ);
        const auto ispin = src.m_basis.ispin(occ);
        if (src.get({!ispin, isite})) {
            DEBUG_ASSERT_EQ(src.site_nocc(isite), 2ul, "should have a double occupation");
            /*
             * this electron on a doubly-occupied site could have been drawn either by this way (uniform branch)
             * or by the preferential branch
             */
            prob = combined_occ_prob(nelec, nelec_doub_occ);
        }
        else {
            /*
             * there is only one way to draw this electron in a singly-occupied site, but the probability must
             * be scaled by the probability that the uniform draw was attempted in the presence of a double
             * occupation
             */
            prob *= 1.0 - m_prob_doub_occ;
        }
        return occ;
    }
    /*
     * preferential branch
     */
    rand = m_prng.draw_uint(nelec_doub_occ);
    const size_t ispin = rand%2;
    const auto isite = doub_occs[rand/2];
    prob = combined_occ_prob(nelec, nelec_doub_occ);
    return src.m_basis.ispinorb(ispin, isite);
}

prob_t exgen::HubbardPreferDoubleOcc::prob_frm(const field::FrmOnv& src, const conn::FrmOnv& conn) const {
    const auto& doub_occs = src.m_decoded.m_doubly_occ_sites.get();
    if (doub_occs.empty()) return prob_uniform(src, conn);
    const auto isite = src.m_basis.isite(conn.m_ann[0]);
    const auto ispin = src.m_basis.ispin(conn.m_ann[0]);
    const auto is_doub_occ = src.get({!ispin, isite});
    if (!is_doub_occ) return prob_uniform(src, conn) * (1-m_prob_doub_occ);
    const auto nelec = src.m_decoded.m_simple_occs.get().size();
    const auto nelec_doub_occ = 2*doub_occs.size();
    const auto nvac = this->valid_adj(isite, src, ispin).size();
    return combined_occ_prob(nelec, nelec_doub_occ)/double(nvac);
}
