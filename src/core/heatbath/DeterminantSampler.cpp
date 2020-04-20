//
// Created by rja on 29/02/2020.
//

#include "DeterminantSampler.h"

DeterminantSampler::DeterminantSampler(const HeatBathSampler &precomputed) :
    m_precomputed(precomputed), m_prng(precomputed.m_prng),
    m_src_det(precomputed.m_h->nsite()), m_anticonn(m_src_det),
    m_occ(m_src_det), m_vac(m_src_det), m_nspinorb(precomputed.m_nbit),
    m_P1(m_nspinorb, 0.0), m_P2_qp(m_nspinorb, 0.0), m_P2_pq(m_nspinorb, 0.0),
    m_P1_aliaser(m_nspinorb, m_prng), m_P2_aliaser(m_nspinorb, m_prng),
    m_P3_aliaser(m_nspinorb, m_prng), m_P4_aliaser(m_nspinorb, m_prng),
    m_single_excitation(m_src_det), m_double_excitation(m_src_det),
    m_single_dst_det(m_src_det), m_double_dst_det(m_src_det) {}


void DeterminantSampler::draw_pq(size_t &p, size_t &q) {

    // (1) draw the first occupied orbital
    p = m_P1_aliaser.draw();
    // (2) iterate over occupied orbital to construct the probability table for the
    // picking of the second occupied orbital. The probability of picking p is set to zero.
    set_P2(m_P2_qp, p);
    m_P2_aliaser.update(m_P2_qp);
    q = m_P2_aliaser.draw();
    ASSERT(q != p);
}

void DeterminantSampler::draw_pq(size_t &p, size_t &q, defs::prob_t &prob) {
    draw_pq(p, q);
    ASSERT(0 < m_P1[p] && m_P1[p] <= 1);
    ASSERT(0 < m_P2_qp[q] && m_P2_qp[q] <= 1);
    prob = m_P1[p] * m_P2_qp[q];
    ASSERT(0 < prob && prob <= 1);
}

void DeterminantSampler::draw_r(const size_t &p, const size_t &q, size_t &r) {
    m_P3_aliaser.update(m_precomputed.m_P3.view(p, q, 0), m_precomputed.m_nbit);
    ASSERT(consts::floats_nearly_equal(m_P3_aliaser.norm(), 1.0, 1e-14));
    r = m_P3_aliaser.draw();
    if (m_src_det.get(r)) r = ~0ul;
}

void DeterminantSampler::draw_r(const size_t &p, const size_t &q, size_t &r, defs::prob_t &prob) {
    draw_r(p, q, r);
    if (r != ~0ul) prob = *m_precomputed.m_P3.view(p, q, r);
    else prob = 0.0;
    ASSERT(0 <= prob && prob <= 1);
}


void DeterminantSampler::draw_pqr(size_t &p, size_t &q, size_t &r) {
    draw_pq(p, q);
    draw_r(p, q, r);
}

void DeterminantSampler::draw_pqr(size_t &p, size_t &q, size_t &r, defs::prob_t &prob) {
    defs::prob_t prob_pq, prob_r;
    draw_pq(p, q, prob_pq);
    draw_r(p, q, r, prob_r);
    prob = prob_pq * prob_r;
    ASSERT(0 <= prob && prob <= 1);
}

void DeterminantSampler::draw_s(const size_t &p, const size_t &q, const size_t &r, size_t &s) {
    // (5) need a double excitation, so choose the second unoccupied
    m_P4_aliaser.update(m_precomputed.m_P4.view(p, q, r, 0), m_nspinorb);
    ASSERT(consts::floats_nearly_equal(m_P4_aliaser.norm(), 1.0, 1e-14));
    s = m_P4_aliaser.draw();
    if (m_src_det.get(s)) s = ~0ul;
}


void DeterminantSampler::draw(size_t &p, size_t &q, size_t &r, size_t &s,
                              defs::prob_t &prob_single, defs::prob_t &prob_double,
                              defs::ham_t &helement_single, defs::ham_t &helement_double) {
    s = ~0ul;
    draw_pqr(p, q, r);
    if (r == ~0ul) {
        prob_single = 0.0;
        prob_double = 0.0;
        return;
    }
    // (4) decide which excitation ranks are to be generated
    m_anticonn.zero();
    m_anticonn.add(p, r);
    m_anticonn.apply(m_src_det);
    helement_single = m_precomputed.m_h->get_element_1(m_anticonn);
    auto htot_rpq = *m_precomputed.m_H_tot.view(p, q, r);
    helement_double = 0.0;

    if (std::abs(helement_single) > htot_rpq) {
        // spawn both single and double
        prob_single = proposal(p, r, helement_single);
        draw_s(p, q, r, s);
        if (s == ~0ul) {
            prob_double = 0.0;
            return;
        }
        helement_double = m_precomputed.m_h->get_element_2(p, q, r, s);
        prob_double = proposal(p, q, r, s, helement_single);
    } else {
        prob_single = std::abs(helement_single) / (htot_rpq + std::abs(helement_single));
        if (m_prng.get().draw_float() < prob_single) {
            // just the single
            prob_double = 0.0;
            prob_single = proposal(p, r, helement_single);
        } else {
            // just the double
            prob_single = 0.0;
            draw_s(p, q, r, s);
            if (s == ~0ul) {
                prob_double = 0.0;
                return;
            }
            helement_double = m_precomputed.m_h->get_element_2(p, q, r, s);
            prob_double = proposal(p, q, r, s, helement_single);
        }
    }
}

void DeterminantSampler::draw() {
    ASSERT(!m_src_det.is_zero()); // call to update method required first
    size_t p, q, r, s;
    defs::prob_t single_prob, double_prob;
    defs::ham_t helement_single, helement_double;
    draw(p, q, r, s, single_prob, double_prob, helement_single, helement_double);

    if (double_prob > 0) {
        m_double_excitation.zero();
        m_double_excitation.add(p, q, r, s);
        m_double_excitation.sort();
        m_double_excitation.apply(m_src_det, m_double_dst_det);
        m_double_prob = double_prob;
        if (single_prob > 0) {
            m_single_excitation.zero();
            m_single_excitation.add(p, r);
            m_single_excitation.apply(m_src_det, m_single_dst_det);
            m_single_prob = single_prob;
            m_outcome = both_excitations;
        } else {
            m_outcome = double_excitation;
        }
    } else {
        if (single_prob > 0) {
            m_single_excitation.zero();
            m_single_excitation.add(p, r);
            m_single_excitation.apply(m_src_det, m_single_dst_det);
            m_single_prob = single_prob;
            m_outcome = single_excitation;
        } else {
            m_outcome = no_excitations;
        }
    }
}

defs::prob_t DeterminantSampler::proposal(const size_t &p, const size_t &r, const defs::ham_t &helement_single) {
    defs::prob_t result = 0.0;
    ASSERT(m_P1[p] > 0.0);
    for (size_t iqp = 0ul; iqp < m_occ.m_nind; ++iqp) {
        const auto &qp = m_occ.m_inds[iqp];
        auto h_tot = *m_precomputed.m_H_tot.view(p, qp, r);
        defs::prob_t p_single = 1.0;
        if (std::abs(helement_single) < h_tot)
            p_single = std::abs(helement_single) / (h_tot + std::abs(helement_single));
        ASSERT(p_single > 0.0);
        result += m_P1[p] * m_P2_qp[qp] *
                  (*m_precomputed.m_P3.view(p, qp, r)) * p_single;
    }
    ASSERT(result > 0.0);
    return result;
}

defs::prob_t DeterminantSampler::proposal(const size_t &p, const size_t &q, const size_t &r, const size_t &s,
                                          const defs::ham_t &helement_single) {
    defs::prob_t result = 0.0;

    auto p_double = [helement_single, this](const size_t &p, const size_t &q, const size_t &r) {
        auto h_tot = *m_precomputed.m_H_tot.view(p, q, r);
        if (std::abs(helement_single) < h_tot) {
            auto tmp = h_tot / (h_tot + std::abs(helement_single));
            ASSERT(tmp > 0);
            return tmp;
        }
        return 1.0;
    };
    ASSERT(*m_precomputed.m_P3.view(p, q, r) >= 0);
    ASSERT(*m_precomputed.m_P3.view(p, q, r) <= 1);
    ASSERT(*m_precomputed.m_P4.view(p, q, r, s) >= 0);
    ASSERT(*m_precomputed.m_P4.view(p, q, r, s) <= 1);

    /*
     * alternative probability given in Appendix B of the reference
     * return 0.25*m_P1[p] * m_P2_qp[q] * p_double(p, q, r) * (*m_precomputed.m_P3.view(p, q, r)) *
     *      (*m_precomputed.m_P4.view(p, q, r, s));
     */

    result += m_P1[p] * m_P2_qp[q] *
              (
                  (*m_precomputed.m_P3.view(p, q, r)) *
                  p_double(p, q, r) *
                  (*m_precomputed.m_P4.view(p, q, r, s))
                  +
                  (*m_precomputed.m_P3.view(p, q, s)) *
                  p_double(p, q, s) *
                  (*m_precomputed.m_P4.view(p, q, s, r))
              );
    set_P2(m_P2_pq, q);

    result += m_P1[q] * m_P2_pq[p] *
              (
                  (*m_precomputed.m_P3.view(q, p, r)) *
                  p_double(q, p, r) *
                  (*m_precomputed.m_P4.view(q, p, r, s))
                  +
                  (*m_precomputed.m_P3.view(q, p, s)) *
                  p_double(q, p, s) *
                  (*m_precomputed.m_P4.view(q, p, s, r))
              );
    ASSERT(result > 0);
    return result;
}

void DeterminantSampler::set_P1(std::vector<defs::prob_t> &P1) {
    P1.assign(P1.size(), 0);
    auto occind = m_occ.m_inds.begin();
    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        if (occind != m_occ.m_inds.end() && p == *occind) {
            P1[p] = *m_precomputed.m_S.view(p);
            occind++;
        }
    }
    prob_utils::normalize(P1);
}


void DeterminantSampler::set_P2(std::vector<defs::prob_t> &P2, const size_t &p) {
    P2.assign(P2.size(), 0);
    auto occind = m_occ.m_inds.begin();
    for (size_t q = 0ul; q < m_nspinorb; ++q) {
        if (occind != m_occ.m_inds.end() && q == *occind) {
            if (q == p) P2[q] = 0.0;
            else P2[q] = *m_precomputed.m_D.view(p, q);
            occind++;
        } else P2[q] = 0.0;
    }
    prob_utils::normalize(P2);
}

void DeterminantSampler::update(const DeterminantElement &det) {
    m_src_det = det;
    m_occ.update(det);
    m_vac.update(det);
    set_P1(m_P1);
    m_P1_aliaser.update(m_P1);
}
