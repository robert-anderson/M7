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
#include "src/sample/Aliaser.h"

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
    NdArray<defs::prob_t, 3> m_P3;
    NdArray<defs::prob_t, 3> m_H_tot;
    NdArray<defs::prob_t, 4> m_P4;

    HeatBathSampler(const Hamiltonian &h);

};


#endif //M7_HEATBATHSAMPLER_H