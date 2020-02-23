//
// Created by Robert John Anderson on 2020-02-22.
//

#include "HeatBathSampler.h"

HeatBathSampler::HeatBathSampler(const AbInitioHamiltonian &h) :
    m_h(h),
    m_nspinorb(h.nspatorb() * 2),
    m_spin_conserving(h.spin_conserving()),
    m_D(m_nspinorb, m_nspinorb),
    m_S(m_nspinorb),
    m_P_tilde_3(m_nspinorb, m_nspinorb, m_nspinorb),
    m_H_tot(m_nspinorb, m_nspinorb, m_nspinorb),
    m_P_tilde_4(m_nspinorb, m_nspinorb, m_nspinorb, m_nspinorb) {

    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        *m_S.view(p) = 0;
        for (size_t q = 0ul; q < m_nspinorb; ++q) {
            if (p == q) continue;
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                for (size_t s = 0ul; s < m_nspinorb; ++s) {
                    *m_D.view(p, q) = std::abs(h.get_element_2(r, s, p, q));
                }
            }
            *m_S.view(p) += *m_D.view(p, q);
        }
    }

    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        for (size_t q = 0ul; q < m_nspinorb; ++q) {
            if (p == q) continue;
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                if (q == r) continue;
                defs::prob_t numerator = 0.0;
                defs::prob_t denominator = 0.0;
                for (size_t sp = 0ul; sp < m_nspinorb; ++sp) {
                    numerator += std::abs(h.int_2().get_phys_antisym(r, sp, p, q));
                    for (size_t rp = 0ul; rp < m_nspinorb; ++rp) {
                        denominator += std::abs(h.int_2().get_phys_antisym(rp, sp, p, q));
                    }
                }
                *m_P_tilde_3.view(p, q, r) = numerator / denominator;
            }
        }
    }

    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        for (size_t q = 0ul; q < m_nspinorb; ++q) {
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                *m_H_tot.view(p, q, r) = 0.0;
                for (size_t s = 0ul; s < m_nspinorb; ++s) {
                    *m_H_tot.view(p, q, r) += std::abs(h.int_2().get_phys_antisym(r, s, p, q));

                    *m_P_tilde_4.view(p, q, r, s) = 0.0;
                    for (size_t sp = 0ul; sp < m_nspinorb; ++sp) {
                        *m_P_tilde_4.view(p, q, r, s) +=
                            std::abs(h.int_2().get_phys_antisym(r, sp, p, q));
                    }
                    *m_P_tilde_4.view(p, q, r, s) =
                        std::abs(h.int_2().get_phys_antisym(r, s, p, q)) /
                        *m_P_tilde_4.view(p, q, r, s);
                }
            }
        }
    }
}

HeatBathSampler::DeterminantSampler HeatBathSampler::sample_excitations(const Determinant &det) {
    return DeterminantSampler(*this, det);
}

HeatBathSampler::DeterminantSampler::DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det) :
    m_precomputed(precomputed), m_det(det),
    m_occinds(det.setinds()), m_noccind(m_occinds.size()),
    m_uncinds(det.clrinds()), m_nuncind(m_uncinds.size()),
    m_P_tilde_1(m_noccind, 0.0),
    m_P_tilde_2(m_noccind, 0.0),
    m_P_tilde_3(det.nspatorb() * 2, 0.0),
    m_P_tilde_4(det.nspatorb() * 2, 0.0),
    m_P_tilde_1_aliaser(m_P_tilde_1) {
    for (size_t ioccind = 0ul; ioccind < m_noccind; ++ioccind) {
        m_P_tilde_1[ioccind] = *m_precomputed.m_S.view(m_occinds[ioccind]);
    }
    prob_utils::normalize(m_P_tilde_1);
}

HeatBathSampler::HeatBathExcitation
HeatBathSampler::DeterminantSampler::draw(PRNG prng) {
    /*
     * This method follows the algorithm detailed in Appendix A of the paper
     */

    // (1) draw the first occupied orbital
    size_t ip = m_occinds[m_P_tilde_1_aliaser.draw(prng)];
    size_t p = m_occinds[ip];

    // (2) iterate over occupied orbital to construct the probability table for the
    // picking of the second occupied orbital. The probability of picking p is set to zero.
    for (size_t iq = 0ul; iq < m_noccind; ++iq) {
        auto q = m_occinds[iq];
        if (iq == ip) m_P_tilde_2[iq] = 0.0;
        else m_P_tilde_2[iq] = *m_precomputed.m_D.view(p, q);
    }
    prob_utils::normalize(m_P_tilde_2);

    Aliaser P_tilde_2_aliaser(m_P_tilde_2);
    size_t iq = P_tilde_2_aliaser.draw(prng);
    size_t q = m_occinds[iq];

    // (3) try to generate the first unoccupied orbital index
    for (size_t r = 0ul; r < m_det.nspatorb() * 2; ++r) {
        if (r == ip || r == iq) m_P_tilde_3[r] = 0.0;
        m_P_tilde_3[r] = *m_precomputed.m_P_tilde_3.view(p, q, r);
    }
    prob_utils::normalize(m_P_tilde_3);
    Aliaser P_tilde_3_aliaser(m_P_tilde_3);
    size_t r = P_tilde_3_aliaser.draw(prng);

    // this draw is made out of all spin orbitals, if the candidate unoccupied orbital
    // we have drawn is in fact occupied, then we have a null excitation
    if (m_det.get(r)) return HeatBathExcitation{Excitation(m_det), Excitation(m_det)};

    // (4) decide which excitation ranks are to be generated
    auto h_rp = m_precomputed.m_h.get_element_1(p, r);
    auto htot_rpq = *m_precomputed.m_H_tot.view(r, p, q);

    defs::prob_t p_single = 0.0;
    defs::prob_t p_double = 0.0;
    if (std::abs(h_rp) > htot_rpq) {
        // spawn both single and double
        p_single = 1.0;
        p_double = 1.0;
    } else {
        p_single = std::abs(h_rp) / (htot_rpq + std::abs(h_rp));
        if (prng.draw_float() < p_single) {
            // just the single
            auto factor = proposal(r, p, p_single);
            return HeatBathExcitation{Excitation(m_det, ip, r, factor), Excitation(m_det)};
        } else {
            p_single = 0.0;
        }
    }

    // (5) need a double excitation, so choose the second unoccupied
    for (size_t s = 0ul; s < m_det.nspatorb() * 2; ++s) {
        if (s == ip || s == iq || s == r) m_P_tilde_3[s] = 0.0;
        m_P_tilde_4[s] = *m_precomputed.m_P_tilde_4.view(p, q, r, s);
    }

    Aliaser P_tilde_4_aliaser(m_P_tilde_4);
    size_t s = P_tilde_4_aliaser.draw(prng);

    auto h_rspq = m_precomputed.m_h.get_element_2(p, q, r, s);

    if (m_det.get(s)) {
        /*
         * no double excitation generated
         */
        if (p_single > 0) {
            return HeatBathExcitation{
                Excitation(m_det, ip, r, proposal(r, p, h_rp)),
                Excitation(m_det)
            };
        } else {
            return HeatBathExcitation{
                Excitation(m_det),
                Excitation(m_det)
            };
        }
    } else {
        if (p_single > 0) {
            return HeatBathExcitation{
                Excitation(m_det, ip, r, proposal(r, p, h_rp)),
                Excitation(m_det, ip, iq, r, s, proposal(ip, iq, r, s, h_rspq))
            };
        } else {
            return HeatBathExcitation{
                Excitation(m_det),
                Excitation(m_det, ip, iq, r, s, proposal(ip, iq, r, s, h_rspq))
            };
        }
    }
}

defs::prob_t HeatBathSampler::DeterminantSampler::proposal(const size_t &ip, const size_t &r, const defs::ham_t &h) {
    defs::prob_t result = 0.0;
    auto p = m_occinds[ip];
    for (auto qp: m_occinds) {
        auto h_tot = *m_precomputed.m_H_tot.view(p, qp, r);
        defs::prob_t p_single = 1.0;
        if (std::abs(h) < h_tot) p_single = std::abs(h) / (h_tot + std::abs(h));
        result += m_P_tilde_1[ip] * m_P_tilde_2[qp] *
                  (*m_precomputed.m_P_tilde_3.view(p, qp, r)) * p_single;
    }
    return result;
}

defs::prob_t
HeatBathSampler::DeterminantSampler::proposal(const size_t &ip, const size_t &iq, const size_t &r, const size_t &s,
                                              const defs::ham_t &h) {
    defs::prob_t result = 0.0;
    auto p = m_occinds[ip];
    auto q = m_occinds[iq];

    auto p_double = [h, this](const size_t &p, const size_t &q, const size_t &r){
        auto h_tot = *m_precomputed.m_H_tot.view(p, q, r);
        if (std::abs(h) < h_tot) return h_tot / (h_tot + std::abs(h));
        return 1.0;
    };

    result += m_P_tilde_1[ip] * m_P_tilde_2[iq] *
              (
                  (*m_precomputed.m_P_tilde_3.view(p, q, r)) *
                  p_double(p, q, r) *
                  (*m_precomputed.m_P_tilde_4.view(p, q, r, s))
                  +
                  (*m_precomputed.m_P_tilde_3.view(p, q, s)) *
                  p_double(p, q, s) *
                  (*m_precomputed.m_P_tilde_4.view(p, q, s, r))
              );
    // repurpose m_P_tilde_2 for the case in which q was picked first
    for (size_t ipp = 0ul; ipp < m_noccind; ++ipp) {
        auto pp = m_occinds[ipp];
        if (ipp == iq) m_P_tilde_2[ipp] = 0.0;
        else m_P_tilde_2[ipp] = *m_precomputed.m_D.view(q, pp);
    }
    prob_utils::normalize(m_P_tilde_2);

    result += m_P_tilde_1[iq] * m_P_tilde_2[ip] *
              (
                  (*m_precomputed.m_P_tilde_3.view(q, p, r)) *
                  p_double(q, p, r) *
                  (*m_precomputed.m_P_tilde_4.view(q, p, r, s))
                  +
                  (*m_precomputed.m_P_tilde_3.view(q, p, s)) *
                  p_double(q, p, s) *
                  (*m_precomputed.m_P_tilde_4.view(q, p, s, r))
              );
    return result;
}
