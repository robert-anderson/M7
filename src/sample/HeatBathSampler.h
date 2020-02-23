//
// Created by Robert John Anderson on 2020-02-22.
//

#ifndef M7_HEATBATHSAMPLER_H
#define M7_HEATBATHSAMPLER_H


#include <src/defs.h>
#include "src/multidim/NdArray.h"
#include "src/integrals/AbInitioHamiltonian.h"
#include "src/enumerators/BitfieldEnumerator.h"
#include "Aliaser.h"

class HeatBathSampler {
public:
    const size_t m_nspinorb;
    const bool m_spin_conserving;
    NdArray<defs::prob_t, 2> m_elec_pair_select_probs;
    NdArray<defs::prob_t, 1> m_elec_select_probs;
    NdArray<defs::prob_t, 3> m_first_hole_select_probs;
    NdArray<defs::prob_t, 3> m_htot;
    NdArray<defs::prob_t, 4> m_double_excitation_probs;

    HeatBathSampler(const AbInitioHamiltonian &h) :
        m_nspinorb(h.nspatorb() * 2),
        m_spin_conserving(h.spin_conserving()),
        m_elec_pair_select_probs(m_nspinorb, m_nspinorb),
        m_elec_select_probs(m_nspinorb),
        m_first_hole_select_probs(m_nspinorb, m_nspinorb, m_nspinorb),
        m_htot(m_nspinorb, m_nspinorb, m_nspinorb),
        m_double_excitation_probs(m_nspinorb, m_nspinorb, m_nspinorb, m_nspinorb) {

        for (size_t p = 0ul; p < m_nspinorb; ++p) {
            *m_elec_select_probs.view(p) = 0;
            for (size_t q = 0ul; q < m_nspinorb; ++q) {
                if (p == q) continue;
                for (size_t r = 0ul; r < m_nspinorb; ++r) {
                    for (size_t s = 0ul; s < m_nspinorb; ++s) {
                        *m_elec_pair_select_probs.view(p, q) = std::abs(h.int_2().get_phys_antisym(r, s, p, q));
                    }
                }
                *m_elec_select_probs.view(p) += *m_elec_pair_select_probs.view(p, q);
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
                    *m_first_hole_select_probs.view(p, q, r) = numerator / denominator;
                }
            }
        }

        for (size_t p = 0ul; p < m_nspinorb; ++p) {
            for (size_t q = 0ul; q < m_nspinorb; ++q) {
                for (size_t r = 0ul; r < m_nspinorb; ++r) {
                    *m_htot.view(p, q, r) = 0.0;
                    for (size_t s = 0ul; s < m_nspinorb; ++s) {
                        *m_htot.view(p, q, r) += std::abs(h.int_2().get_phys_antisym(r, s, p, q));

                        *m_double_excitation_probs.view(p, q, r, s) = 0.0;
                        for (size_t sp = 0ul; sp < m_nspinorb; ++sp) {
                            *m_double_excitation_probs.view(p, q, r, s) +=
                                std::abs(h.int_2().get_phys_antisym(r, sp, p, q));
                        }
                        *m_double_excitation_probs.view(p, q, r, s) =
                            std::abs(h.int_2().get_phys_antisym(r, s, p, q)) /
                            *m_double_excitation_probs.view(p, q, r, s);
                    }
                }
            }
        }
    }

    class DeterminantSampler{
        const HeatBathSampler &m_precomputed;
        const Determinant &m_det;
        const defs::inds m_occinds;
        const size_t m_noccind;
        const defs::inds m_uncinds;
        const size_t m_nuncind;
        std::vector<defs::prob_t> m_first_elec_probs;
        const Aliaser m_first_elec_aliaser;
        std::vector<defs::prob_t> m_second_elec_probs;
        std::vector<defs::prob_t> m_first_hole_probs;
        std::vector<defs::prob_t> m_second_hole_probs;
        DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det):
        m_precomputed(precomputed), m_det(det),
        m_occinds(det.setinds()), m_noccind(m_occinds.size()),
        m_uncinds(det.clrinds()), m_nuncind(m_uncinds.size()),
        m_first_elec_probs(m_noccind, 0.0),
        m_second_elec_probs(m_noccind-1, 0.0),
        m_first_hole_probs(m_nuncind, 0.0),
        m_second_hole_probs(m_nuncind-1, 0.0),
        m_first_elec_aliaser(m_first_elec_probs)
        {
            size_t ioccind = 0ul;
            for (auto occind : m_occinds){
                m_first_elec_probs[ioccind] = *m_precomputed.m_elec_select_probs.view(m_occinds[ioccind]);
                ioccind++;
            }
        }

        void draw(PRNG prng, size_t &elec_1, size_t &elec_2, size_t &hole_1, size_t &hole_2){
            elec_1 = m_occinds[m_first_elec_aliaser.draw(prng)];
            size_t ioccind = 0ul;
            for (auto occind : m_occinds){
                if (occind==elec_1) continue;
                m_second_elec_probs[ioccind] = *m_precomputed.m_elec_pair_select_probs.view(elec_1, occind);
                ioccind++;
            }
            Aliaser second_elec_aliaser(m_second_elec_probs);
            elec_2 = m_occinds[second_elec_aliaser.draw(prng)];

            size_t iuncind = 0ul;
            for (auto uncind : m_uncinds){
                if (m_precomputed.m_spin_conserving && m_det.orbspin(elec_1)!=m_det.orbspin(uncind)) continue;
                //m_first_hole_probs[iuncind] = *m_precomputed.m
                iuncind++;
            }

        }
    };

    NdArray<defs::prob_t, 1> get_first_elec_probs(const Determinant &det){
        auto setbits = DeterminantSetEnumerator(det).enumerate();
        auto nsetbit = setbits.size();
    }

};


#endif //M7_HEATBATHSAMPLER_H
