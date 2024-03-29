/**
 * @file DenseTcHamiltonian.cpp
 * @author jph
 * @brief test file for exact diagonalisation of dense transcorrelated fermion Hamiltonians.
 * @date 2022-05-30
 *
 */

#include <test_core/defs.h>
#include "M7_lib/linalg/DenseHamiltonian.h"
#include "M7_lib/io/Symlink.h"
#include "M7_lib/hamiltonian/frm/TcFrmHam.h"
#include "M7_lib/util/Sort.h"

#ifdef ENABLE_TCHINT

/**
 * @brief exact diagonalisation of dense TC Be atom in 6-31G basis
 *
 */
TEST(DenseTcFermionHamiltonian, TcBe631G) {
    AssetSymlink tcdump("TC_Be_6-31G/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_6-31G/FCIDUMP", "FCIDUMP");
    TcFrmHam ham_src({"FCIDUMP"});
    Hamiltonian gham(&ham_src);
    DenseHamiltonian ham(gham);
    // non-symmetric real matrices have (in general) complex eigenvalues
    v_t<std::complex<ham_comp_t>> evals;
    dense::diag(ham, evals); // dense diagonalisation
    // unlike hermitian case, we need to sort the eigenvalue array by the real part (lowest first)
    sort::inplace(evals, false, false);
    // std::cout.precision(17);
    // std::cout << evals[0] << std::endl; // print ground state
    // compare to NECI calculation
    ASSERT_NEAR_EQ(evals[0].real(), -14.666331930789127);
    // check that all eigenvalues have zero imaginary component
    for (unsigned int i=0; i < evals.size(); ++i) {
        ASSERT_NEAR_ZERO(evals[i].imag());
    }
}

#endif // ENABLE_TCHINT
