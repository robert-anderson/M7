//
// Created by RJA on 20/11/2020.
//

#include "UniformSingles.h"

UniformSingles::UniformSingles(const FermionHamiltonian *ham, PRNG &prng) :
        FermionExcitationGenerator(ham, prng, 1) {}

bool UniformSingles::draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv, const OccupiedOrbitals &occ,
                           const VacantOrbitals &vac, defs::prob_t &prob, defs::ham_t &helem,
                           conn::AsFermionOnv &anticonn) {
    size_t i, a, ia;
    size_t ncases;
    if (m_spin_conserving) {
        size_t nalpha = src_fonv.nalpha();
        size_t nbeta = m_nelec - nalpha;
        size_t nalpha_cases = nalpha * (m_h->nsite() - nalpha);
        size_t nbeta_cases = nbeta * (m_h->nsite() - nbeta);
        ncases = nalpha_cases + nbeta_cases;
        ia = m_prng.draw_uint(ncases);
        if (ia < nalpha_cases) {
            integer_utils::inv_rectmap(i, a, m_h->nsite() - nalpha, ia);
            ASSERT(i < nalpha);
            ASSERT(a < m_h->nsite() - nalpha);
        } else {
            ia -= nalpha_cases;
            integer_utils::inv_rectmap(i, a, m_h->nsite() - nbeta, ia);
            // skip the occupied alphas
            i += nalpha;
            // skip the unoccupied alphas
            a += m_h->nsite() - nalpha;
            ASSERT(i >= nalpha);
            ASSERT(i < m_nelec);
            ASSERT(a >= m_h->nsite() - nalpha);
            ASSERT(a < 2 * m_h->nsite() - m_nelec);
        }
    } else {
        ncases = m_nelec * (2 * m_h->nsite() - m_nelec);
        ia = m_prng.draw_uint(ncases);
        integer_utils::inv_rectmap(i, a, 2 * m_h->nsite() - m_nelec, ia);
    }
    i = occ.m_inds[i];
    a = vac.m_inds[a];
#ifndef NDEBUG
    if (m_spin_conserving) {
        if (i < m_h->nsite()) ASSERT(a < m_h->nsite())
        else ASSERT(a >= m_h->nsite())
    }
#endif
    anticonn.zero();
    anticonn.add(i, a);
    anticonn.apply(src_fonv, dst_fonv);
    prob = 1.0 / (defs::prob_t) (ncases);
    helem = m_h->get_element_1(anticonn);
    return !consts::float_nearly_zero(helem, 1e-12);
}
