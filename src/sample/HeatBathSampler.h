//
// Created by Robert John Anderson on 2020-02-22.
//

#ifndef M7_HEATBATHSAMPLER_H
#define M7_HEATBATHSAMPLER_H


#include <src/defs.h>
#include <src/fermion/Excitation.h>
#include "src/multidim/NdArray.h"
#include "src/integrals/AbInitioHamiltonian.h"
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
    const AbInitioHamiltonian &m_h;
    const size_t m_nspinorb;
    const bool m_spin_conserving;
    NdArray<defs::prob_t, 2> m_D;
    NdArray<defs::prob_t, 1> m_S;
    NdArray<defs::prob_t, 3> m_P_tilde_3;
    NdArray<defs::prob_t, 3> m_H_tot;
    NdArray<defs::prob_t, 4> m_P_tilde_4;

    HeatBathSampler(const AbInitioHamiltonian &h);

    struct HeatBathExcitation {
        Excitation m_single;
        Excitation m_double;
    };

    class DeterminantSampler{
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

        HeatBathExcitation draw(PRNG prng);

    private:
        defs::prob_t proposal(const size_t &p, const size_t &r, const defs::ham_t &h){
            defs::prob_t out = 0.0;
            for (auto qp: m_occinds){
                auto htot = m_precomputed.m_H_tot.view(p, qp, r);
                defs::prob_t p_single = 1.0;
                if ()
                out+=
                    m_P_tilde_1[p]*
                    m_P_tilde_2[qp]*
                    (*m_precomputed.m_P_tilde_3.view(p, qp, r))*
                    prob
            }
        }
    };

    DeterminantSampler sample_excitations(const Determinant &det){
        return DeterminantSampler(*this, det);
    }

};


#endif //M7_HEATBATHSAMPLER_H
