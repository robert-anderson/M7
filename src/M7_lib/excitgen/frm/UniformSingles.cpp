//
// Created by anderson on 04/04/2022.
//

#include "UniformSingles.h"

bool UniformSingles::draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig_utils::ex_single, "this excitation generator is only suitable for exsig 1100");
    auto spin_conserving = m_h.m_kramers_attrs.m_conserving_singles;
    if (spin_conserving) return draw_spin_conserve_fn(m_prng, src, prob, conn);
    return draw_spin_nonconserve_fn(m_prng, src, prob, conn);
}

size_t UniformSingles::approx_nconn() const {
    const auto& hs = m_h.m_hs;
    if (hs.ms2_conserved()) {
        return hs.m_nelec_alpha * hs.m_nvac_alpha + hs.m_nelec_beta * hs.m_nvac_beta;
    } else {
        return hs.m_nelec * hs.m_nvac;
    }
}

bool UniformSingles::draw_spin_conserve_fn(PRNG &prng, const field::FrmOnv &src,
                                           defs::prob_t &prob, conn::FrmOnv &conn) {
    const auto &nonempty_pairs = src.m_decoded.m_nonempty_pair_labels.get();
    if (nonempty_pairs.empty()) return false;
    const auto label = nonempty_pairs[prng.draw_uint(nonempty_pairs.size())];
    const auto &occs = src.m_decoded.m_spin_sym_occs.get()[label];
    const auto &vacs = src.m_decoded.m_spin_sym_vacs.get()[label];

    const auto nocc = occs.size();
    const auto nvac = vacs.size();
    DEBUG_ASSERT_TRUE(nocc, "should have emitted null excitation if number of occupied orbitals is zero");
    DEBUG_ASSERT_TRUE(nvac, "should have emitted null excitation if number of vacant orbitals is zero");

    size_t ia = prng.draw_uint(nocc * nvac);
    size_t i, a;
    integer_utils::inv_rectmap(i, a, nvac, ia);
    DEBUG_ASSERT_LT(i, nocc, "drawn occupied OOB");
    DEBUG_ASSERT_LT(a, nvac, "drawn vacant OOB");

    conn.m_ann.set(occs[i]);
    conn.m_cre.set(vacs[a]);
    DEBUG_ASSERT_TRUE(conn.kramers_conserve(), "spin not conserved");
    prob = 1.0 / (nonempty_pairs.size() * nocc * nvac);
    return true;
}

bool UniformSingles::draw_spin_nonconserve_fn(PRNG &prng, const field::FrmOnv &src,
                                              defs::prob_t &prob, conn::FrmOnv &conn) {
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto &vacs = src.m_decoded.m_simple_vacs.get();
    const auto nelec = occs.size();
    DEBUG_ASSERT_EQ(nelec, src.m_hs.m_nelec, "unexpected number of electrons");
    const auto ncases = nelec * src.m_hs.m_nvac;
    auto ia = prng.draw_uint(ncases);
    size_t i, a;
    integer_utils::inv_rectmap(i, a, src.m_sites.m_nspinorb - nelec, ia);
    conn.m_ann.set(occs[i]);
    conn.m_cre.set(vacs[a]);
    return true;
}
