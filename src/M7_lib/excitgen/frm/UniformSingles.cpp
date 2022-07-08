//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "UniformSingles.h"

bool exgen::UniformSingles::draw_frm(uint_t exsig, const field::FrmOnv& src, prob_t& prob, conn::FrmOnv& conn) {
    DEBUG_ONLY(exsig);
    DEBUG_ASSERT_EQ(exsig, exsig::ex_single, "this excitation generator is only suitable for exsig 1100");
    auto spin_conserving = m_h.m_kramers_attrs.m_conserving_singles;
    if (spin_conserving) return draw_spin_conserve_fn(m_prng, src, prob, conn);
    return draw_spin_nonconserve_fn(m_prng, src, prob, conn);
}

bool exgen::UniformSingles::draw_spin_conserve_fn(PRNG& prng, const field::FrmOnv& src,
                                           prob_t& prob, conn::FrmOnv& conn) {
    const auto& nonempty_pairs = src.m_decoded.m_nonempty_pair_labels.get();
    if (nonempty_pairs.empty()) return false;
    const auto label = nonempty_pairs[prng.draw_uint(nonempty_pairs.size())];
    const auto& occs = src.m_decoded.m_spin_sym_occs.get()[label];
    const auto& vacs = src.m_decoded.m_spin_sym_vacs.get()[label];

    const auto nocc = occs.size();
    const auto nvac = vacs.size();
    DEBUG_ASSERT_TRUE(nocc, "should have emitted null excitation if number of occupied orbitals is zero");
    DEBUG_ASSERT_TRUE(nvac, "should have emitted null excitation if number of vacant orbitals is zero");

    uint_t ia = prng.draw_uint(nocc * nvac);
    uint_t i, a;
    integer::inv_rectmap(i, a, nvac, ia);
    DEBUG_ASSERT_LT(i, nocc, "drawn occupied OOB");
    DEBUG_ASSERT_LT(a, nvac, "drawn vacant OOB");

    conn.m_ann.set(occs[i]);
    conn.m_cre.set(vacs[a]);
    DEBUG_ASSERT_TRUE(conn.kramers_conserve(), "spin not conserved");
    prob = 1.0 / (nonempty_pairs.size() * nocc * nvac);
    return true;
}

bool exgen::UniformSingles::draw_spin_nonconserve_fn(PRNG& prng, const field::FrmOnv& src,
                                              prob_t& prob, conn::FrmOnv& conn) {
    const auto& occs = src.m_decoded.m_simple_occs.get();
    const auto& vacs = src.m_decoded.m_simple_vacs.get();
    const auto nelec = occs.size();
    const auto nvac = src.m_basis.m_nspinorb - nelec;
    const auto ncases = nelec * nvac;
    auto ia = prng.draw_uint(ncases);
    uint_t i, a;
    integer::inv_rectmap(i, a, nvac, ia);
    conn.m_ann.set(occs[i]);
    conn.m_cre.set(vacs[a]);
    prob = 1.0 / ncases;
    return true;
}

prob_t exgen::UniformSingles::prob_spin_conserve_fn(const field::FrmOnv& src, const conn::FrmOnv& conn) {
    const auto& nonempty_pairs = src.m_decoded.m_nonempty_pair_labels.get();
    if (nonempty_pairs.empty()) return 0.0;
    const auto label = src.m_decoded.m_spin_sym_occs.label(conn.m_cre[0]);
    const auto& occs = src.m_decoded.m_spin_sym_occs.get()[label];
    const auto& vacs = src.m_decoded.m_spin_sym_vacs.get()[label];
    return 1.0 / (nonempty_pairs.size() * occs.size() * vacs.size());
}

prob_t exgen::UniformSingles::prob_spin_nonconserve_fn(const field::FrmOnv& src, const conn::FrmOnv&) {
    const auto nocc = src.m_decoded.m_simple_occs.get().size();
    const auto nvac = src.m_basis.m_nspinorb - nocc;
    return 1.0/(nocc*nvac);
}

prob_t exgen::UniformSingles::prob_fn(const field::FrmOnv& src, const conn::FrmOnv& conn, bool spin_conserve) {
    return spin_conserve ? prob_spin_conserve_fn(src, conn) : prob_spin_nonconserve_fn(src, conn);
}

prob_t exgen::UniformSingles::prob_frm(const field::FrmOnv& src, const conn::FrmOnv& conn) const {
    auto spin_conserving = m_h.m_kramers_attrs.m_conserving_singles;
    if (spin_conserving) return prob_spin_conserve_fn(src, conn);
    return prob_spin_nonconserve_fn(src, conn);
}

uint_t exgen::UniformSingles::approx_nconn(uint_t, sys::Particles particles) const {
    const auto& elecs = particles.m_frm;
    sys::frm::Sector sector(m_h.m_basis, particles.m_frm);
    if (elecs.m_ms2.conserve()){
        return elecs.m_nalpha * sector.m_nvac_alpha + elecs.m_nbeta * sector.m_nvac_beta;
    } else {
        return elecs * sector.m_nvac;
    }
}