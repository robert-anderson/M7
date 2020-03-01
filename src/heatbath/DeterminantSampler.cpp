//
// Created by rja on 29/02/2020.
//

#include "DeterminantSampler.h"


DeterminantSampler::DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det) :
    m_precomputed(precomputed), m_det(det),
    m_occinds(det.setinds()), m_noccind(m_occinds.size()),
    m_uncinds(det.clrinds()), m_nuncind(m_uncinds.size()),
    m_P1(make_P1(precomputed, m_occinds, m_noccind)),
    m_P2_qp(det.nspatorb() * 2, 0.0),
    m_P2_pq(det.nspatorb() * 2, 0.0),
    m_P1_aliaser(m_P1) {}

void DeterminantSampler::draw_pq(PRNG &prng, size_t &p, size_t &q) {

    // (1) draw the first occupied orbital
    p = m_P1_aliaser.draw(prng);
    // (2) iterate over occupied orbital to construct the probability table for the
    // picking of the second occupied orbital. The probability of picking p is set to zero.
    make_P2(m_P2_qp, p);
    Aliaser P_tilde_2_aliaser(m_P2_qp);
    q = P_tilde_2_aliaser.draw(prng);
    assert(q != p);
}

void DeterminantSampler::draw_pq(PRNG &prng, size_t &p, size_t &q, defs::prob_t &prob) {
    draw_pq(prng, p, q);
    assert(0<m_P1[p] && m_P1[p]<=1);
    assert(0<m_P2_qp[q] && m_P2_qp[q]<=1);
    prob = m_P1[p] * m_P2_qp[q];
    assert(0<prob && prob<=1);
}

void DeterminantSampler::draw_r(PRNG &prng, const size_t &p, const size_t &q, size_t &r){
    Aliaser aliaser(m_precomputed.m_P3.view(p, q, 0), m_precomputed.m_nspinorb);
    assert(consts::floats_nearly_equal(aliaser.norm(), 1.0, 1e-14));
    r = aliaser.draw(prng);
    if (m_det.get(r)) r = ~0ul;
}

void DeterminantSampler::draw_r(PRNG &prng, const size_t &p, const size_t &q, size_t &r, defs::prob_t &prob){
    draw_r(prng, p, q, r);
    if (r!=0ul) prob = *m_precomputed.m_P3.view(p, q, r);
    else prob = 0.0;
    assert(0<=prob && prob<=1);
}


void DeterminantSampler::draw_pqr(PRNG &prng, size_t &p, size_t &q, size_t &r){
    draw_pq(prng, p, q);
    draw_r(prng, p, q, r);
}

void DeterminantSampler::draw_pqr(PRNG &prng, size_t &p, size_t &q, size_t &r, defs::prob_t &prob){
    defs::prob_t prob_pq, prob_r;
    draw_pq(prng, p, q, prob_pq);
    draw_r(prng, p, q, r, prob_r);
    prob = prob_pq*prob_r;
    assert(0<=prob && prob<=1);
}

void DeterminantSampler::draw_s(PRNG &prng, const size_t &p, const size_t &q, const size_t &r, size_t &s) {
    // (5) need a double excitation, so choose the second unoccupied
    Aliaser aliaser(m_precomputed.m_P4.view(p, q, r, 0), m_precomputed.m_nspinorb);
    assert(consts::floats_nearly_equal(aliaser.norm(), 1.0, 1e-14));
    s = aliaser.draw(prng);
    if (m_det.get(s)) s = ~0ul;
}


void DeterminantSampler::draw(PRNG &prng, size_t &p, size_t &q, size_t &r, size_t &s,
                                   defs::prob_t &prob_single, defs::prob_t &prob_double,
                                   defs::ham_t &helement_single, defs::ham_t &helement_double){
    s = ~0ul;
    draw_pqr(prng, p, q, r);
    if (r==~0ul){
        prob_single = 0.0;
        prob_double = 0.0;
        return;
    }
    // (4) decide which excitation ranks are to be generated
    helement_single = m_precomputed.m_h.get_element_1(m_det, p, r);
    auto htot_rpq = *m_precomputed.m_H_tot.view(p, q, r);
    helement_double = 0.0;

    if (std::abs(helement_single) > htot_rpq) {
        // spawn both single and double
        prob_single = proposal(p, r, helement_single);
        draw_s(prng, p, q, r, s);
        if (s==~0ul) {
            prob_double = 0.0; return;
        }
        helement_double = m_precomputed.m_h.get_element_2(p, q, r, s);
        prob_double = proposal(p, q, r, s, helement_double);
    } else {
        prob_single = std::abs(helement_single) / (htot_rpq + std::abs(helement_single));
        if (prng.draw_float() < prob_single) {
            // just the single
            prob_double = 0.0;
            prob_single = proposal(r, p, helement_single);
        } else {
            // just the double
            prob_single = 0.0;
            draw_s(prng, p, q, r, s);
            if (s==~0ul) {
                prob_double = 0.0; return;
            }
            helement_double = m_precomputed.m_h.get_element_2(p, q, r, s);
            prob_double = proposal(p, q, r, s, helement_double);
        }
    }
}

HeatBathExcitation DeterminantSampler::draw(PRNG &prng) {
    size_t p, q, r, s;
    defs::prob_t prob_single, prob_double;
    defs::ham_t helement_single, helement_double;
    draw(prng, p, q, r, s, prob_single, prob_double, helement_single, helement_double);


    if (r != ~0ul && s != ~0ul)
        return HeatBathExcitation{
            Excitation(m_det),
            Excitation(m_det, p, q, r, s,
                       m_P1[p] * m_P2_qp[q] *(*m_precomputed.m_P3.view(p, q, r)) *
                       (*m_precomputed.m_P4.view(p, q, r, s)),
                helement_double)
        };
    else {
        return HeatBathExcitation{
            Excitation(m_det),
            Excitation(m_det)
        };
    }

    if (prob_double>0) {
        if (prob_single > 0)
            return HeatBathExcitation{
                Excitation(m_det, p, q, prob_single, helement_single),
                Excitation(m_det, p, q, r, s, prob_double, helement_double)
            };
        else
            return HeatBathExcitation{
                Excitation(m_det),
                Excitation(m_det, p, q, r, s, prob_double, helement_double)
            };
    } else {
        if (prob_single > 0)
            return HeatBathExcitation{
                Excitation(m_det, p, q, prob_single, helement_single),
                Excitation(m_det)
            };
        else
            return HeatBathExcitation{
                Excitation(m_det),
                Excitation(m_det)
            };
    }
}

defs::prob_t DeterminantSampler::proposal(const size_t &p, const size_t &r, const defs::ham_t &h) {
    defs::prob_t result = 0.0;
    for (auto qp: m_occinds) {
        auto h_tot = *m_precomputed.m_H_tot.view(p, qp, r);
        defs::prob_t p_single = 1.0;
        if (std::abs(h) < h_tot) p_single = std::abs(h) / (h_tot + std::abs(h));
        result += m_P1[p] * m_P2_qp[qp] *
                  (*m_precomputed.m_P3.view(p, qp, r)) * p_single;
    }
    return result;
}

defs::prob_t DeterminantSampler::proposal(const size_t &p, const size_t &q, const size_t &r, const size_t &s,
                                          const defs::ham_t &h) {
    defs::prob_t result = 0.0;

/*
    auto p_double = [h, this](const size_t &p, const size_t &q, const size_t &r) {
        auto h_tot = *m_precomputed.m_H_tot.view(p, q, r);
        if (std::abs(h) < h_tot) return h_tot / (h_tot + std::abs(h));
        return 1.0;
    };
    */
    auto p_double = [h, this](const size_t &p, const size_t &q, const size_t &r) {
        return 1.0;

        auto h_tot = *m_precomputed.m_H_tot.view(p, q, r);
        if (std::abs(h) < h_tot) {
            auto tmp = h_tot / (h_tot + std::abs(h));
            assert(tmp > 0);
            return tmp;
        }
        return 1.0;
    };
    assert(*m_precomputed.m_P3.view(p, q, r)>=0);
    assert(*m_precomputed.m_P3.view(p, q, r)<=1);
    assert(*m_precomputed.m_P4.view(p, q, r, s)>=0);
    assert(*m_precomputed.m_P4.view(p, q, r, s)<=1);
    //return m_P1[p] * m_P2_qp[q] * p_double(p, q, r) * (*m_precomputed.m_P3.view(p, q, r)) *
    //       (*m_precomputed.m_P4.view(p, q, r, s));
    return m_P1[p] * m_P2_qp[q] *(*m_precomputed.m_P3.view(p, q, r)) *
           (*m_precomputed.m_P4.view(p, q, r, s));



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
    make_P2(m_P2_pq, q);

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
    assert(result > 0);
    return result;
}

std::vector<defs::prob_t> DeterminantSampler::make_P1(const HeatBathSampler &precomputed, const defs::inds &occinds,
                                                      const size_t &noccind) {
    const auto n = precomputed.m_S.nelement();
    std::vector<defs::prob_t> result(n, 0);
    auto occind = occinds.begin();
    assert(occinds.begin()+noccind==occinds.end());
    for (size_t p = 0ul; p < n; ++p) {
        if (occind!=occinds.end() && p == *occind) {
            result[p] = *precomputed.m_S.view(p);
            occind++;
        }
    }
    prob_utils::normalize(result);
    return result;
}


void DeterminantSampler::make_P2(std::vector<defs::prob_t> &P2, const size_t &p) {
    const auto n = P2.size();
    auto occind = m_occinds.begin();
    for (size_t q = 0ul; q < n; ++q) {
        if (occind!=m_occinds.end() && q == *occind) {
            if (q == p) P2[q] = 0.0;
            else P2[q] = *m_precomputed.m_D.view(p, q);
            occind++;
        } else P2[q] = 0.0;
    }
    prob_utils::normalize(P2);
}