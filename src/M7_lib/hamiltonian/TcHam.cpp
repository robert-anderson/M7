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
 * @return ham_t the corresponding L matrix element
 */
ham_t TcHam::get_lmat_coeff(uint_t a, uint_t b, uint_t c, uint_t i,
                                  uint_t j, uint_t k) const {
#ifdef ENABLE_TCHINT
    const int ia = a+1, ib = b+1, ic = c+1, ii = i+1, ij = j+1, ik = k+1;  // convert to int
    return three_body_matel(&ia, &ib, &ic, &ii, &ij, &ik);
    // return 2.71;
#else
    dtype::unused(a, b, c, i, j, k);
    return 0.0;
#endif
}
