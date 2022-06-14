/**
 * @file DenseTcHamiltonian.cpp
 * @author jph
 * @brief test file for exact diagonalisation of dense transcorrelated fermion Hamiltonians.
 * @date 2022-05-30
 *
 */

#include <gtest/gtest.h>
#include "M7_lib/linalg/DenseHamiltonian.h"
#include "M7_lib/io/Symlink.h"
#include "M7_lib/hamiltonian/TcFrmHam.h"

#ifdef ENABLE_TCHINT

/**
 * @brief exact diagonalisation of dense TC Be atom in 6-31G basis
 *
 */
 /*
TEST(DenseTcFermionHamiltonian, TcBe631G) {
    AssetSymlink tcdump("TC_Be_6-31G/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_6-31G/FCIDUMP", "FCIDUMP");
    TcFrmHam ham_src("FCIDUMP", false, 0);
    // DenseHamiltonian ham(ham_src);
    // std::vector<defs::ham_t> evals;
    // dense::diag(ham, evals);
    // TODO stub
    ASSERT(false);
}
  */

#endif // ENABLE_TCHINT