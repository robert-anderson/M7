//
// Created by rja on 21/07/22.
//

#include "HolsteinUniform.h"
#include "M7_lib/hamiltonian/frmbos/HolsteinLadderHam.h"

exgen::HolsteinUniform::HolsteinUniform(const FrmBosHam& h, PRNG& prng, OpSig exsig) :
        FrmBosExcitGen(h, prng, {exsig}, "uniform Holstein "+str_t(exsig.nbos_cre() ? "creation" : "annhilation")){
    REQUIRE_TRUE(h.is<HolsteinLadderHam>(), "holstein excit gen requires holstein hamiltonian");
    REQUIRE_EQ(h.m_basis.m_frm.m_nsite, h.m_basis.m_bos.m_nmode,
               "holstein excit gen assumes one boson mode per fermion site");
}

uint_t exgen::HolsteinUniform::approx_nconn(OpSig /*exsig*/, sys::Particles particles) const {
    return particles.m_frm;
}

bool exgen::HolsteinUniform0010::draw_frmbos(
        OpSig exsig, const field::FrmBosOnv& src, prob_t& prob, conn::FrmBosOnv& conn) {
    DEBUG_ONLY(exsig);
    DEBUG_ASSERT_EQ(exsig, opsig::c_0010, "invalid exsig specified");
    const auto& occ_cutoff = m_h.m_basis.m_bos.m_occ_cutoff;
    if (!m_h.m_basis.m_bos.m_occ_cutoff) return false;
    const auto &occs = src.m_frm.m_decoded.m_simple_occs.get();

    auto imode = m_h.m_basis.m_frm.isite(occs[m_prng.draw_uint(occs.size())]);
    // attempt to generate an ONV with an additional boson occupying the mode at imode
    size_t curr_occ = src.m_bos[imode];
    DEBUG_ASSERT_LE(curr_occ, occ_cutoff, "current occupation of selected mode exceeds cutoff");

    // doubly-occupied sites are twice as likely to be drawn
    prob = src.m_frm.site_nocc(imode);
    prob /= occs.size();

    if (curr_occ == occ_cutoff) return false;

    conn.clear();
    conn.m_bos.m_cre.set(imode);
    return true;
}

prob_t exgen::HolsteinUniform0010::prob_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn) const {
    const auto imode = conn.m_bos.m_cre[0].m_imode;
    const auto nelec = src.m_frm.m_decoded.m_simple_occs.get().size();
    return src.m_frm.site_nocc(imode)/double(nelec);
}

bool exgen::HolsteinUniform0001::draw_frmbos(
        OpSig exsig, const field::FrmBosOnv& src, prob_t& prob, conn::FrmBosOnv& conn) {
    DEBUG_ONLY(exsig)
    DEBUG_ASSERT_EQ(exsig, opsig::c_0001, "invalid exsig specified");
    if (!m_h.m_basis.m_bos.m_occ_cutoff) return false;

    const auto &occs = src.m_decoded.m_occ_sites_nonzero_bosons.get();
    if (occs.empty()) return false;

    auto imode = occs[m_prng.draw_uint(occs.size())];
    // attempt to generate a "de-excited" ONV with one less boson occupying the mode at imode
    DEBUG_ASSERT_LE(src.m_bos[imode], m_h.m_basis.m_bos.m_occ_cutoff,
                    "current occupation of selected mode exceeds cutoff");
    DEBUG_ASSERT_TRUE(src.m_bos[imode], "selected boson mode should have non-zero occupation");

    prob = 1.0/occs.size();

    conn.clear();
    conn.m_bos.m_ann.set(imode);
    return true;
}

prob_t exgen::HolsteinUniform0001::prob_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& /*conn*/) const {
    return 1.0/double(src.m_decoded.m_occ_sites_nonzero_bosons.get().size());
}