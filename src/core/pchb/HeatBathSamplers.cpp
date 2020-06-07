//
// Created by rja on 09/05/2020.
//

#include "HeatBathSamplers.h"

HeatBathSamplers::HeatBathSamplers(const Hamiltonian *h, PrivateStore<PRNG> &prng):
    ExcitationGenerator(h, prng){
    std::vector<defs::prob_t> weights(m_norb_pair, 0.0);
    size_t ij = 0ul;
    size_t ab = 0ul;
    for (size_t i = 0ul; i < m_norb; ++i) {
        for (size_t j = 0ul; j < i; ++j) {
            weights.assign(m_norb_pair, 0.0);
            ab = 0ul;
            for (size_t a = 0ul; a < m_norb; ++a) {
                for (size_t b = 0ul; b < a; ++b) {
                    //if (a!=i && a!=j && b!=i && b!=j) { !TODO why does this restriction fail?
                        weights[ab] = std::abs(m_h->get_element_2(i, j, a, b));
                    //}
                    ++ab;
                }
            }
            m_pick_ab_given_ij.emplace_back(weights, prng);
            ++ij;
        }
    }
    ASSERT(ij == m_norb_pair);
    ASSERT(ab == m_norb_pair);
}

bool HeatBathSamplers::draw_single(const DeterminantElement &src_det, DeterminantElement &dst_det,
                                   const OccupiedOrbitals &occ, const VacantOrbitals &vac, defs::prob_t &prob,
                                   defs::ham_t &helem, AntisymConnection &anticonn) {
    size_t i, a, ia;
    size_t ncases;
    if (m_spin_conserving) {
        size_t nalpha = src_det.nalpha();
        size_t nbeta = m_nelec - nalpha;
        size_t nalpha_cases = nalpha * (m_h->nsite() - nalpha);
        size_t nbeta_cases = nbeta * (m_h->nsite() - nbeta);
        ncases = nalpha_cases + nbeta_cases;
        ia = m_prng.get().draw_uint(ncases);
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
            ASSERT(i < m_nelec);
            ASSERT(a >= m_h->nsite() - nalpha);
            ASSERT(a < 2*m_h->nsite()-m_nelec);
        }
    } else {
        ncases = m_nelec * (2 * m_h->nsite() - m_nelec);
        ia = m_prng.get().draw_uint(ncases);
        integer_utils::inv_rectmap(i, a, 2 * m_h->nsite() - m_nelec, ia);
    }
    i = occ.m_inds[i];
    a = vac.m_inds[a];
#ifndef NDEBUG
    if (m_spin_conserving) {
        if (i<m_h->nsite()) ASSERT(a<m_h->nsite())
        else ASSERT(a>=m_h->nsite())
    }
#endif
    anticonn.zero();
    anticonn.add(i, a);
    anticonn.apply(src_det, dst_det);
    prob = 1.0 / (defs::prob_t) (ncases);
    helem = m_h->get_element_1(anticonn);
    return !consts::float_nearly_zero(helem, 1e-12);
}

bool HeatBathSamplers::draw_double(const DeterminantElement &src_det, DeterminantElement &dst_det,
                                   const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem,
                                   AntisymConnection &anticonn) {
    // just draw uniform ij TODO! int weighted ij
    // return false if invalid excitation generated, true otherwise
    size_t i, j, a, b;
    size_t ij = m_prng.get().draw_uint(m_nelec_pair);
    integer_utils::inv_strigmap(j, i, ij);
    // i and j are positions in the occ list, convert to orb inds:
    i = occ.m_inds[i];
    j = occ.m_inds[j];
    ASSERT(std::any_of(occ.m_inds.cbegin(), occ.m_inds.cbegin()+occ.m_nind, [&i](const size_t &k) {return k == i;}));
    ASSERT(std::any_of(occ.m_inds.cbegin(), occ.m_inds.cbegin()+occ.m_nind, [&j](const size_t &k) {return k == j;}));
    ASSERT(i<j);

    ij = integer_utils::strigmap(j, i); // i and j are orbital indices
    size_t ab = m_pick_ab_given_ij[ij].draw();
    integer_utils::inv_strigmap(b, a, ab); // a and b are orbital indices
    //ASSERT(i!=a && i!=b && j!=a && j!=b)

    auto either_vac_in_array = [&a, &b](const size_t &k) {return k == a || k == b;};

    if (std::any_of(occ.m_inds.cbegin(), occ.m_inds.cbegin()+occ.m_nind, either_vac_in_array)) {
        return 0;
    }

    anticonn.zero();
    anticonn.add(i, j, a, b);
    anticonn.apply(src_det, dst_det);
    helem = m_h->get_element_2(anticonn);
    prob = std::abs(helem) / (m_pick_ab_given_ij[ij].norm() * m_nelec_pair);
    if (consts::float_nearly_zero(prob, 1e-14)) {
        return false;
    }
    return true;
}
