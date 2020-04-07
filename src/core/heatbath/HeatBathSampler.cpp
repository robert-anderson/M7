//
// Created by Robert John Anderson on 2020-02-22.
//

#include "HeatBathSampler.h"

HeatBathSampler::HeatBathSampler(const Hamiltonian &h) :
        m_h(h),
        m_nbit(h.nsite() * 2),
        m_spin_conserving(h.spin_conserving()),
        m_D(m_nbit, m_nbit),
        m_S(m_nbit),
        m_P3(m_nbit, m_nbit, m_nbit),
        m_H_tot(m_nbit, m_nbit, m_nbit),
        m_P4(m_nbit, m_nbit, m_nbit, m_nbit) {

    for (size_t p = 0ul; p < m_nbit; ++p) {
        *m_S.view(p) = 0;
        for (size_t q = 0ul; q < m_nbit; ++q) {
            if (p == q) continue;
            for (size_t r = 0ul; r < m_nbit; ++r) {
                for (size_t s = 0ul; s < m_nbit; ++s) {
                    *m_D.view(p, q) += std::abs(h.get_element_2(p, q, r, s));
                }
            }
            *m_S.view(p) += *m_D.view(p, q);
        }
    }

    for (size_t p = 0ul; p < m_nbit; ++p) {
        for (size_t q = 0ul; q < m_nbit; ++q) {
            if (p == q) continue;
            defs::ham_comp_t denom = 0.0;
            for (size_t r = 0ul; r < m_nbit; ++r) {
                if (p == r || q == r) continue;
                for (size_t s = 0ul; s < m_nbit; ++s) {
                    if (p == s || q == s || r == s) continue;
                    *m_P3.view(p, q, r) += std::abs(h.get_element_2(p, q, r, s));
                }
                denom += *m_P3.view(p, q, r);
            }
            for (size_t r = 0ul; r < m_nbit; ++r) {
                if (!consts::float_is_zero(denom)) *m_P3.view(p, q, r) /= denom;
            }
            auto norm = std::accumulate(m_P3.view(p, q, 0), m_P3.view(p, q, 0) + m_nbit, 0.0);
            assert(consts::float_is_zero(norm) || consts::floats_nearly_equal(norm, 1.0));
        }
    }

    for (size_t p = 0ul; p < m_nbit; ++p) {
        for (size_t q = 0ul; q < m_nbit; ++q) {
            if (p == q) continue;
            for (size_t r = 0ul; r < m_nbit; ++r) {
                if (p == r || q == r) continue;
                *m_H_tot.view(p, q, r) = 0.0;
                for (size_t s = 0ul; s < m_nbit; ++s) {
                    if (p == s || q == s || r == s) continue;
                    *m_H_tot.view(p, q, r) += std::abs(h.get_element_2(p, q, r, s));
                }
                for (size_t s = 0ul; s < m_nbit; ++s) {
                    if (p == s || q == s || r == s) continue;
                    if (!consts::float_is_zero(*m_H_tot.view(p, q, r)))
                        *m_P4.view(p, q, r, s) =
                                std::abs(h.get_element_2(p, q, r, s)) / (*m_H_tot.view(p, q, r));
                }
                auto norm = std::accumulate(m_P4.view(p, q, r, 0),
                                            m_P4.view(p, q, r, 0) + m_nbit, 0.0);
                assert(consts::float_is_zero(norm) || consts::floats_nearly_equal(norm, 1.0));
            }
        }
    }
}
