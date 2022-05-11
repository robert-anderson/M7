/**
 * @file TcHam.cpp
 * @author jph
 * @brief Implementation file for general transcorrelated Hamiltonian
 * @date 2022-05-11
 *
 */

#include "TcHam.h"

TcHam::TcHam() {
#ifdef ENABLE_TCHINT
    tchint_init();
#else
    MPI_ABORT("Compiled without TCHInt support.");
#endif
}

TcHam::~TcHam() {
#ifdef ENABLE_TCHINT
    tchint_finalize();
#else
    MPI_ABORT("Compiled without TCHInt support.");
#endif
}

/**
 * @brief wrapper to TCHInt's L matrix (6-index integrals)
 *
 * @return defs::ham_t the corresponding L matrix element
 */
defs::ham_t TcHam::get_lmat_coeff(size_t a, size_t b, size_t c, size_t i,
                                  size_t j, size_t k) const {
#ifdef ENABLE_TCHINT
    const int ia = i, ib = b, ic = c, ii = i, ij = j, ik = k;  // convert to int
    return three_body_matel(&ia, &ib, &ic, &ii, &ij, &ik);
#else
    return 0.0
#endif
}
