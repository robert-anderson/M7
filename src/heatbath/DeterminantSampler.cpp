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
        m_P3(det.nspatorb() * 2, 0.0),
        m_P4(det.nspatorb() * 2, 0.0),
        m_P1_aliaser(m_P1) {}

HeatBathExcitation DeterminantSampler::draw(PRNG &prng) {
    /*
     * This method follows the algorithm detailed in Appendix A of the paper
     */

    // (1) draw the first occupied orbital
    size_t p = m_P1_aliaser.draw(prng);

    // (2) iterate over occupied orbital to construct the probability table for the
    // picking of the second occupied orbital. The probability of picking p is set to zero.
    make_P2(m_P2_qp, p);

    Aliaser P_tilde_2_aliaser(m_P2_qp);
    size_t q = P_tilde_2_aliaser.draw(prng);
    assert(q != p);

    // (3) try to generate the first unoccupied orbital index
    make_P3(m_P3, p, q);
    Aliaser P_tilde_3_aliaser(m_P3);
    size_t r = P_tilde_3_aliaser.draw(prng);

    // this draw is made out of all spin orbitals, if the candidate unoccupied orbital
    // we have drawn is in fact occupied, then we have a null excitation
    if (m_det.get(r)) return HeatBathExcitation{Excitation(m_det), Excitation(m_det)};

    // (4) decide which excitation ranks are to be generated
    auto h_rp = m_precomputed.m_h.get_element_1(m_det, p, r);
    auto htot_rpq = *m_precomputed.m_H_tot.view(p, q, r);

    defs::prob_t p_single = 0.0;
    /*

if (std::abs(h_rp) > htot_rpq) {
    // spawn both single and double
    p_single = 1.0;
    p_double = 1.0;
} else {
    p_single = std::abs(h_rp) / (htot_rpq + std::abs(h_rp));
    if (prng.draw_float() < p_single) {
        // just the single
        auto prob = proposal(r, p, h_rp);
        return HeatBathExcitation{Excitation(m_det, p, r, h_rp, prob), Excitation(m_det)};
    } else {
        p_single = 0.0;
    }
}
 */

    // (5) need a double excitation, so choose the second unoccupied
    make_P4(m_P4, p, q, r);

    Aliaser P_tilde_4_aliaser(m_P4);
    size_t s = P_tilde_4_aliaser.draw(prng);

    auto h_rspq = m_precomputed.m_h.get_element_2(p, q, r, s);
    return HeatBathExcitation{
            Excitation(m_det),
            Excitation(m_det, p, q, r, s, h_rspq, proposal(p, q, r, s, h_rspq))
    };



    if (m_det.get(s)) {
        /*
         * no double excitation generated
         */

        if (p_single > 0) {
            /*
            assert(0);
            return HeatBathExcitation{
                    Excitation(m_det, p, r, h_rp, proposal(r, p, h_rp)),
                    Excitation(m_det)
            };*/
        } else {
            assert(0);
            return HeatBathExcitation{
                    Excitation(m_det),
                    Excitation(m_det)
            };
        }
    } else {
        if (p_single > 0) {
            assert(0);
            return HeatBathExcitation{
                    Excitation(m_det, p, r, h_rp, proposal(r, p, h_rp)),
                    Excitation(m_det, p, q, r, s, h_rspq, proposal(p, q, r, s, h_rspq))
            };
        } else {
            return HeatBathExcitation{
                    Excitation(m_det),
                    Excitation(m_det, p, q, r, s, h_rspq, proposal(p, q, r, s, h_rspq))
            };
        }
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
    return m_P1[p] * m_P2_qp[q] * p_double(p, q, r) * (*m_precomputed.m_P3.view(p, q, r)) *
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
    for (size_t p = 0ul; p < n; ++p) {
        if (p == *occind) {
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
        if (q == *occind) {
            if (q == p) P2[q] = 0.0;
            else P2[q] = *m_precomputed.m_D.view(p, q);
            occind++;
        } else P2[q] = 0.0;
    }
    prob_utils::normalize(P2);
}


void DeterminantSampler::make_P3(std::vector<defs::prob_t> &P3, const size_t &p, const size_t &q) {
    const auto n = P3.size();
    for (size_t r = 0ul; r < n; ++r) {
        if (r == p || r == q) m_P3[r] = 0.0;
        else P3[r] = *m_precomputed.m_P3.view(p, q, r);
    }
    prob_utils::normalize(P3);
}

void DeterminantSampler::make_P4(std::vector<defs::prob_t> &P4, const size_t &p, const size_t &q, const size_t &r) {
    const auto n = P4.size();
    for (size_t s = 0ul; s < n; ++s) {
        if (s == p || s == q || s==r) m_P4[s] = 0.0;
        else P4[s] = *m_precomputed.m_P4.view(p, q, r, s);
    }
    prob_utils::normalize(P4);
}