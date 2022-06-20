/**
 * @file TranscorrelatedBosonHamiltonian.cpp
 * @author jph
 * @brief test file for transcorrelated Boson Hamiltonians
 * @note
 * Coming up with benchmarks for transcorrelated Bosons is hairy business, so
 * at the moment the benchmarks are simply from carefully combining matrix
 * elements manually for a small system
 *
 * @par
 * We are using FCIDUMP and TCDUMP files generated for an ab initio fermion
 * Hamiltonian, but treat it as if it were bosonic
 *
 * @date 2022-06-01
 *
 */

#include <M7_lib/hamiltonian/bos/TcBosHam.h>
#include <M7_lib/io/Symlink.h>
