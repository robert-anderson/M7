//
// Created by anderson on 04/04/2022.
//

#include "UniformSingles2.h"

bool UniformSingles2::draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                                 conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig_utils::ex_single, "this excitation generator is only suitable for exsig 1100");

    auto spin_conserving = m_h.m_kramers_attrs.m_conserving_singles;
    size_t i, a, ia;
    size_t ncases;
    if (spin_conserving) {
        const auto &nonempty_pairs = src.m_decoded.m_nonempty_pair_labels.get();
        if (nonempty_pairs.empty()) return false;
        const auto label = nonempty_pairs[m_prng.draw_uint(nonempty_pairs.size())];
        const auto &occs = src.m_decoded.m_spin_sym_occs.get()[label];
        const auto &vacs = src.m_decoded.m_spin_sym_vacs.get()[label];

        const auto nocc = occs.size();
        const auto nvac = vacs.size();
        DEBUG_ASSERT_TRUE(nocc, "should have emitted null excitation if number of occupied orbitals is zero");
        DEBUG_ASSERT_TRUE(nvac, "should have emitted null excitation if number of vacant orbitals is zero");

        ncases = nonempty_pairs.size() * nocc * nvac;
        ia = m_prng.draw_uint(nocc * nvac);
        integer_utils::inv_rectmap(i, a, nvac, ia);
        DEBUG_ASSERT_LT(i, nocc, "drawn occupied case OOB");
        DEBUG_ASSERT_LT(a, nvac, "drawn vacant case OOB");
        i = occs[i];
        a = vacs[a];
    } else {
        const auto nspinorb = 2 * m_h.m_nsite;
        ncases = m_h.m_nelec * (nspinorb - m_h.m_nelec);
        ia = m_prng.draw_uint(ncases);
        integer_utils::inv_rectmap(i, a, nspinorb - m_h.m_nelec, ia);
        i = src.m_decoded.m_simple_occs.get()[i];
        a = src.m_decoded.m_simple_vacs.get()[a];
    }
#ifndef NDEBUG
    if (spin_conserving) {
        if (i < m_h.m_nsite) { DEBUG_ASSERT_LT(a, m_h.m_nsite, "spin not conserved"); }
        else { DEBUG_ASSERT_GE(a, m_h.m_nsite, "spin not conserved"); }
    }
#endif
    conn.set(i, a);
    prob = 1.0 / ncases;
    return true;
}

size_t UniformSingles2::approx_nconn() const {
    auto spin_conserving = m_h.m_kramers_attrs.m_conserving_singles;
    if (spin_conserving) {
        return 2 * (m_h.m_nelec / 2) * ((2*m_h.m_nsite - m_h.m_nelec) / 2);
    } else {
        return m_h.m_nelec * (2*m_h.m_nsite - m_h.m_nelec);
    }
}
