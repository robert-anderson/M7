/**
 * @file DenseTcBosonHamiltonian.cpp
 * @author jph
 * @brief test file for exact diagonalisation of dense transcorrelated Boson Hamiltonian
 *
 * @par
 * at the moment, we have no benchmarks for this so the tests are all stubs
 * @date 2022-06-15
 *
 */

#include <test_core/defs.h>
#include "M7_lib/linalg/DenseHamiltonian.h"
#include "M7_lib/io/Symlink.h"
#include "M7_lib/hamiltonian/bos/TcBosHam.h"
#include "M7_lib/util/Sort.h"

#ifdef ENABLE_TCHINT

/**
 * @brief exact diagonalisation of dense TC "bosonic Be atom" in 6-31G basis
 *
 */
TEST(DenseTcBosonHamiltonian, TcBe631G) {
    AssetSymlink tcdump("TC_Be_6-31G/TCDUMP", "TCDUMP");
    AssetSymlink bosdump("TC_Be_6-31G/BOSDUMP", "BOSDUMP");
    AssetSymlink fcidump("TC_Be_6-31G/FCIDUMP", "FCIDUMP");
    TcBosHam ham_src({"BOSDUMP"}, sys::bos::c_max_occ);
    Hamiltonian gham(&ham_src);
    DenseHamiltonian ham(gham);
    // non-symmetric real matrices have (in general) complex eigenvalues
    v_t<std::complex<ham_comp_t>> evals;
    dense::diag(ham, evals); // dense diagonalisation
    // unlike hermitian case, we need to sort the eigenvalue array by the real part (lowest first)
    sort::inplace(evals, true, false);
    std::cout.precision(17);
    std::cout << evals << std::endl; // print ground state
    // compare to NECI calculation
    ASSERT(false);
}

#endif // ENABLE_TCHINT
