/**
 * @file TranscorrelatedHamiltonian.cpp
 * @author jph
 * @brief test file for general transcorrelated Hamiltonians
 * @date 2022-05-24
 *
 */

#include <M7_lib/hamiltonian/TcHam.h> // what's being tested
#include <gtest/gtest.h>

#ifdef ENABLE_TCHINT
/**
 * @brief check values of some get_coeff_element
 *
 */
TEST(TranscorrelatedHamiltonian, get_coeff3300_test) {
    // General Hamiltonian to be tested
    TcHam ham;
    // remember 1-based: becomes (1,1,1,1,1,1)
    ASSERT_FLOAT_EQ(ham.get_lmat_coeff(0,0,0,0,0,0)/3., -0.13564812863025142E-001);
    // (2, 2, 2, 2, 4, 4)
    ASSERT_FLOAT_EQ(ham.get_lmat_coeff(1,1,3,1,3,1)/3., -0.74599783654206124E-004 );
    // 1 2 1 1 2 5 -> 1 1 2 1 5 2 -> 2 2 3 2 6 3
    ASSERT_FLOAT_EQ(ham.get_lmat_coeff(1,2,1,5,2,1)/3., 0.80828847497582549E-003);
    // 3 2 6 2 7 3 -> 2 3 6 7 2 3 -> 2 2 3 7 3 6 (=0 in this case)
    ASSERT_FLOAT_EQ(ham.get_lmat_coeff(2,1,5,1,6,2)/3., 0.0);
    ASSERT_FLOAT_EQ(ham.get_lmat_coeff(1,1,2,3,4,5)/3., 0.0);
    // 2 3 4 5 6 7
    ASSERT_FLOAT_EQ(ham.get_lmat_coeff(1,2,3,4,5,6)/3., 0.40443217097041021E-009);
}

#endif // ENABLE_TCHINT
