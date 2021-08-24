//
// Created by RJA on 20/11/2020.
//

#include "UniformSingles.h"

UniformSingles::UniformSingles(const Hamiltonian &ham, PRNG &prng) :
        FrmExcitGen(ham, prng, 1) {}

bool UniformSingles::draw(const size_t& exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                          defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) {
    size_t i, a, ia;
    size_t ncases;
    if (m_spin_conserving) {
        const auto& nonempty_pairs = orbs.nonempty_pair_labels(src);
        if (nonempty_pairs.empty()) return false;
        const auto label = nonempty_pairs[m_prng.draw_uint(nonempty_pairs.size())];
        const auto& occs = orbs.occ(src)[label];
        const auto& vacs = orbs.vac(src)[label];

        const auto nocc = occs.size();
        const auto nvac = vacs.size();
        DEBUG_ASSERT_TRUE(nocc, "should have emitted null excitation if number of occupied orbitals is zero");
        DEBUG_ASSERT_TRUE(nvac, "should have emitted null excitation if number of vacant orbitals is zero");

        ncases = nonempty_pairs.size()*nocc*nvac;
        ia = m_prng.draw_uint(nocc*nvac);
        integer_utils::inv_rectmap(i, a, nvac, ia);
        DEBUG_ASSERT_LT(i, nocc, "drawn occupied case OOB");
        DEBUG_ASSERT_LT(a, nvac, "drawn vacant case OOB");
        i = occs[i];
        a = occs[a];
    } else {
        ncases = m_nelec * (2 * m_h.nsite() - m_nelec);
        ia = m_prng.draw_uint(ncases);
        integer_utils::inv_rectmap(i, a, 2 * m_h.nsite() - m_nelec, ia);
        i = orbs.occ(src).m_flat[i];
        a = orbs.vac(src).m_flat[a];
    }
#ifndef NDEBUG
    if (m_spin_conserving) {
        if (i < m_h.nsite()) ASSERT(a < m_h.nsite())
        else ASSERT(a >= m_h.nsite())
    }
#endif
    conn.add(i, a);
    prob = 1.0 / ncases;
    helem = m_h.m_frm.get_element_1100(src, conn);
    return !consts::float_nearly_zero(helem, 1e-12);
}

size_t UniformSingles::approx_nconn() const {
    if(m_spin_conserving) {
        return 2*(m_nelec/2)*((m_nspinorb-m_nelec)/2);
    } else {
        return m_nelec*(m_nspinorb-m_nelec);
    }
}
