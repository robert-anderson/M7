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
#include <M7_lib/util/utils.h>

struct TcHam {
    // Would there be any members? I guess not really
    TcHam();  // constructor

    ~TcHam();  // destructor

    // protected:
    defs::ham_t get_lmat_coeff(size_t a, size_t b, size_t c, size_t i, size_t j,
                               size_t k) const;
};

#endif  // M7_TCHAM_H
