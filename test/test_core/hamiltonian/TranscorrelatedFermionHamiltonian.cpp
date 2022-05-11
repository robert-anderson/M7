/**
 * @file TranscorrelatedFermionHamiltonian.cpp
 * @author jph
 * @brief test file for transcorrelated Fermion Hamiltonians
 * @date 2022-05-03
 *
 */

#include <gtest/gtest.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/field/Mbf.h>
#include "M7_lib/caches/DecodedDeterminants.h"

// TODO:
// [ ] put appropriate TCDUMP(s) into the assets folder
// [ ] get_coeff2200 check for non-Hermiticity (should be handled fine)
// [ ] get_coeff3300 same as get_element3300 up to sign
//          test with and without the parity
// [ ] contracted elements (get_element{00,11,22}00)
//
// [ ] maybe also do a "ui test" as done in TCHINT

#ifdef ENABLE_TCHINT
/**
 * @brief checks if get_coeff_element3300 and get_coeff_element3300 are the same
 *        up to parity
 *
 */
TEST(TranscorrelatedFermionHamiltonian, coeff_element3300_parity) {
    // TODO stub
}
#endif // ENABLE_TCHINT
