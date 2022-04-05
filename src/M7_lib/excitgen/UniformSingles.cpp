//
// Created by RJA on 20/11/2020.
//

#include "UniformSingles.h"

UniformSingles::UniformSingles(const Hamiltonian &h, PRNG &prng) :
        FrmExcitGen(h, prng, exsig_utils::ex_single) {}

bool UniformSingles::draw_frm(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                          defs::prob_t &prob, conn::FrmOnv &conn) {
    size_t i, a, ia;
    size_t ncases;
    if (m_spin_conserving) {
        const auto &nonempty_pairs = orbs.nonempty_pair_labels(src);
        if (nonempty_pairs.empty()) return false;
        const auto label = nonempty_pairs[m_prng.draw_uint(nonempty_pairs.size())];
        const auto &occs = orbs.occ(src)[label];
        const auto &vacs = orbs.vac(src)[label];

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
        ncases = m_nelec * (m_h.m_bd.m_nspinorb - m_nelec);
        ia = m_prng.draw_uint(ncases);
        integer_utils::inv_rectmap(i, a, m_h.m_bd.m_nspinorb - m_nelec, ia);
        i = orbs.occ(src).m_flat[i];
        a = orbs.vac(src).m_flat[a];
    }
#ifndef NDEBUG
    if (m_spin_conserving) {
        if (i < m_h.m_bd.m_nsite) { DEBUG_ASSERT_LT(a, m_h.m_bd.m_nsite, "spin not conserved"); }
        else { DEBUG_ASSERT_GE(a, m_h.m_bd.m_nsite, "spin not conserved"); }
    }
#endif
    conn.m_ann.set(i);
    conn.m_cre.set(a);
    prob = 1.0 / ncases;
    return true;
}

size_t UniformSingles::approx_nconn() const {
    if (m_spin_conserving) {
        return 2 * (m_nelec / 2) * ((m_bd.m_nspinorb - m_nelec) / 2);
    } else {
        return m_nelec * (m_bd.m_nspinorb - m_nelec);
    }
}

std::string UniformSingles::description() const {
    return "uniform";
}
