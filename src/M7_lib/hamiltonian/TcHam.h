/**
 * @file TcHam.h
 * @author jph
 * @brief general transcorrelated Hamiltonian
 *  This is for multiple inheritance (i.e. this does not have any declarations
 *  coeffs or matels)
 * @date 2022-05-09
 *
 */

#ifndef M7_TCHAM_H
#define M7_TCHAM_H

#ifdef ENABLE_TCHINT
#include "tchint.h"
#endif

#include <M7_lib/parallel/MPIAssert.h>

struct TcHam {
    // Would there be any members? I guess not really
    TcHam();  // constructor

    ~TcHam();  // destructor

    // protected:
    ham_t get_lmat_coeff(uint_t a, uint_t b, uint_t c, uint_t i, uint_t j, uint_t k) const;
};

#endif  // M7_TCHAM_H
