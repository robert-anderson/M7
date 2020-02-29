//
// Created by Robert John Anderson on 2020-02-22.
//

#ifndef M7_HEATBATHSAMPLER_H
#define M7_HEATBATHSAMPLER_H


#include <src/defs.h>
#include <src/fermion/Excitation.h>
#include "src/multidim/NdArray.h"
#include "src/hamiltonian/Hamiltonian.h"
#include "src/enumerators/BitfieldEnumerator.h"
#include "Aliaser.h"

/*
 * An implementation of the Heat Bath Excitation generator of A. A. Holmes et al
 *      J. Chem. Theory Comput. 2016, 12, 1561−1571
 *
 * Notational conventions are adopted from the paper.
 *
 * HeatBathSampler is generated in shared memory at the outset of the calculation,
 * and remains constant.
 *
 * Instances of DeterminantSampler are generated in private memory for each occupied
 * determinant.
 *
 * DeterminantSampler's draw method yields a HeatBathExcitation instance, which either
 * holds a single Excitation, a double Excitation, neither, or both. The proper use of
 * this method is to call it a number of time roughly proportional to the weight of the
 * determinant in question.
 */


class HeatBathSampler {
public:
    const Hamiltonian &m_h;
    const size_t m_nspinorb;
    const bool m_spin_conserving;
    NdArray<defs::prob_t, 2> m_D;
    NdArray<defs::prob_t, 1> m_S;
    NdArray<defs::prob_t, 3> m_P_tilde_3;
    NdArray<defs::prob_t, 3> m_H_tot;
    NdArray<defs::prob_t, 4> m_P_tilde_4;

    HeatBathSampler(const Hamiltonian &h);

    struct HeatBathExcitation {
        Excitation m_single;
        Excitation m_double;
    };

    class DeterminantSampler {
        const HeatBathSampler &m_precomputed;
        const Determinant &m_det;
        const defs::inds m_occinds;
        const size_t m_noccind;
        const defs::inds m_uncinds;
        const size_t m_nuncind;
        std::vector<defs::prob_t> m_P_tilde_1;
        const Aliaser m_P_tilde_1_aliaser;
        std::vector<defs::prob_t> m_P_tilde_2;
        std::vector<defs::prob_t> m_P_tilde_3;
        std::vector<defs::prob_t> m_P_tilde_4;

    public:
        DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det);

        std::vector<defs::prob_t> make_m_P_tilde_1(const HeatBathSampler &precomputed,
                                    const defs::inds &occinds, const size_t &noccind) {
            std::vector<defs::prob_t> result(noccind, 0ul);
            for (size_t ioccind = 0ul; ioccind < noccind; ++ioccind) {
                result[ioccind] = *precomputed.m_S.view(occinds[ioccind]);
            }
            prob_utils::normalize(result);
            return result;
        }

        HeatBathExcitation draw(PRNG &prng);

    private:
        defs::prob_t proposal(const size_t &ip, const size_t &r, const defs::ham_t &h);

        defs::prob_t proposal(const size_t &ip, const size_t &iq,
                              const size_t &r, const size_t &s,
                              const defs::ham_t &h);
    };

    DeterminantSampler sample_excitations(const Determinant &det) const;

};


#endif //M7_HEATBATHSAMPLER_H
