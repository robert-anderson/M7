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
#include "M7_lib/hamiltonian/frm/TcFrmHam.h"

#ifdef ENABLE_TCHINT

/**
 * @brief exact diagonalisation of dense TC Be atom in 6-31G basis
 *
 */
TEST(DenseTcFermionHamiltonian, TcBe631G) {
    AssetSymlink tcdump("TC_Be_6-31G/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_6-31G/FCIDUMP", "FCIDUMP");
    TcFrmHam ham_src({"FCIDUMP"}, false);
    Hamiltonian gham(&ham_src);
    DenseHamiltonian ham(gham);
    // non-symmetric real matrices have complex eigenvalues
    std::vector<std::complex<defs::ham_t>> evals;
    dense::diag(ham, evals); // dense diagonalisation
    // unlike hermitian case, we need to sort the eigenvalue array by the real part (lowest first)
    sort_utils::inplace(evals, false, false);
    std::cout << evals[0] << std::endl; // print ground state
    // compare to NECI calculation
    ASSERT_DOUBLE_EQ(evals[0].real(), -14.6663);
    ASSERT_DOUBLE_EQ(evals[0].imag(), 0.0);
}

#endif // ENABLE_TCHINT
