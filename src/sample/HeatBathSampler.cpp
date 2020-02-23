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
                    *m_D.view(p, q) = std::abs(h.int_2().get_phys_antisym(r, s, p, q));
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

HeatBathSampler::DeterminantSampler::DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det) :
    m_precomputed(precomputed), m_det(det),
    m_occinds(det.setinds()), m_noccind(m_occinds.size()),
    m_uncinds(det.clrinds()), m_nuncind(m_uncinds.size()),
    m_P_tilde_1(m_noccind, 0.0),
    m_P_tilde_2(m_noccind - 1, 0.0),
    m_P_tilde_3(det.nspatorb()*2, 0.0),
    m_P_tilde_4(det.nspatorb()*2, 0.0),
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
        m_P_tilde_3[r] = *m_precomputed.m_P_tilde_3.view(p, q, r);
    }
    prob_utils::normalize(m_P_tilde_3);
    Aliaser P_tilde_3_aliaser(m_P_tilde_3);
    size_t ir = P_tilde_3_aliaser.draw(prng);
    size_t r = m_occinds[ir];

    // this draw is made out of all spin orbitals, if the candidate unoccupied orbital
    // we have drawn is in fact occupied, then we have a null excitation
    if (m_det.get(r)) return HeatBathExcitation{Excitation(m_det), Excitation(m_det)};

    // (4) decide which excitation ranks are to be generated
    auto h_rp = m_precomputed.m_h.get_element_1(p, r);
    auto htot_rpq = *m_precomputed.m_H_tot.view(r, p, q);

    defs::prob_t p_single = 0.0;
    defs::prob_t p_double = 0.0;
    defs::ham_t factor;
    if (std::abs(h_rp) > htot_rpq) {
        // spawn both single and double
        p_single = 1.0;
        p_double = 1.0;
    }
    else {
        p_single = std::abs(h_rp) / (htot_rpq + std::abs(h_rp));
        if (prng.draw_float() < p_single) {
            // just the single
            factor = proposal(r, p, p_single);
            return HeatBathExcitation{Excitation(m_det, p, r, factor), Excitation(m_det)};
        }
        else {
            // just the double
        }
    }

}
