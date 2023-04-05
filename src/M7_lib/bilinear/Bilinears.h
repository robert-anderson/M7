//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_BILINEARS_H
#define M7_BILINEARS_H

#include "M7_lib/connection/OpSig.h"
/**
 * Projection onto a trial wavefunction is sufficient for the estimation of many-body expectation values if the operator
 * in question commutes with the Hamiltonian. When the operator does not commute with the Hamiltonian, the projection
 * by a trial wavefunction is just an approximation called a "mixed estimator".
 *
 * In this general case, the only way to obtain estimates which are exact in the infinite sampling limit is to use the
 * wavefunction itself as the trial wavefunction, and thus we obtain estimates of a quantity which is "bilinear" in the
 * wavefunction.
 *
 * This is statistically problematic when the WF is stochastically propagated, since products of correlated walker
 * populations introduce a systematic bias. This is overcome by replicating walker populations for each root.
 *
 * Two different types of multidimensional bilinear estimator are currently implemented.
 *  1. reduced density matrices (RDMs)
 *  2. spectral moments (SpecMoms)
 *
 * In M7, RDMs are defined as the expectation value of a normal ordered product of second quantized (SQ) operators and
 * may include fermionic spin orbital indices or bosonic modes as free indices.
 *
 * SpecMoms are defined as matrix elements of some power of the Hamiltonian between perturbed wavefunctions, where the
 * perturbing operators have arbitrary rank but are limited to fermionic spin orbital indices.
 */

namespace bilinears {
    /**
     * @param string
     *  excitation signature or fermion operator rank defined in the configuration document
     * @return
     *  integer excitation signature
     */
    OpSig parse_exsig(const str_t &string);
    /**
     * @param strings
     *  excitation signatures or fermion operator ranks defined in the configuration document
     * @return
     *  integer excitation signatures
     */
    v_t<OpSig> parse_exsigs(const strv_t &strings);

    bool in_parsed_exsigs(const OpSig& opsig, const strv_t &strings);

}

#endif //M7_BILINEARS_H