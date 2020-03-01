//
// Created by Robert John Anderson on 2020-02-22.
//

#include "HeatBathSampler.h"

HeatBathSampler::HeatBathSampler(const Hamiltonian &h) :
    m_h(h),
    m_nspinorb(h.nsite() * 2),
    m_spin_conserving(h.spin_conserving()),
    m_D(m_nspinorb, m_nspinorb),
    m_S(m_nspinorb),
    m_P3(m_nspinorb, m_nspinorb, m_nspinorb),
    m_H_tot(m_nspinorb, m_nspinorb, m_nspinorb),
    m_P4(m_nspinorb, m_nspinorb, m_nspinorb, m_nspinorb) {

    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        *m_S.view(p) = 0;
        for (size_t q = 0ul; q < m_nspinorb; ++q) {
            if (p == q) continue;
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                for (size_t s = 0ul; s < m_nspinorb; ++s) {
                    *m_D.view(p, q) += std::abs(h.get_element_2(p, q, r, s));
                }
            }
            *m_S.view(p) += *m_D.view(p, q);
        }
    }

    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        for (size_t q = 0ul; q < m_nspinorb; ++q) {
            if (p == q) continue;
            defs::ham_comp_t denom = 0.0;
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                if (p == r || q == r) continue;
                for (size_t s = 0ul; s < m_nspinorb; ++s) {
                    if (p == s || q == s || r == s) continue;
                    *m_P3.view(p, q, r) += std::abs(h.get_element_2(p, q, r, s));
                }
                denom += *m_P3.view(p, q, r);
            }
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                *m_P3.view(p, q, r) /= denom;
            }
            assert(consts::floats_nearly_equal(
                std::accumulate(m_P3.view(p, q, 0), m_P3.view(p, q, 0) + m_nspinorb, 0.0), 1.0));
        }
    }

    for (size_t p = 0ul; p < m_nspinorb; ++p) {
        for (size_t q = 0ul; q < m_nspinorb; ++q) {
            if (p == q) continue;
            for (size_t r = 0ul; r < m_nspinorb; ++r) {
                if (p == r || q == r) continue;
                *m_H_tot.view(p, q, r) = 0.0;
                for (size_t s = 0ul; s < m_nspinorb; ++s) {
                    if (p == s || q == s || r == s) continue;
                    *m_H_tot.view(p, q, r) += std::abs(h.get_element_2(p, q, r, s));
                }
                for (size_t s = 0ul; s < m_nspinorb; ++s) {
                    if (p == s || q == s || r == s) continue;
                    *m_P4.view(p, q, r, s) =
                        std::abs(h.get_element_2(p, q, r, s)) / (*m_H_tot.view(p, q, r));
                }
                assert(consts::floats_nearly_equal(
                    std::accumulate(m_P4.view(p, q, r, 0),
                                    m_P4.view(p, q, r, 0) + m_nspinorb, 0.0), 1.0));
            }
        }
    }
}