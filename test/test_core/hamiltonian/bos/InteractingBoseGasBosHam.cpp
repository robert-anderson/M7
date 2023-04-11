//
// Created by Robert J. Anderson on 08/03/2022.
//

#include "test_core/defs.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"
#include "M7_lib/hamiltonian/Hamiltonian.h"

TEST(InteractingBosGasBosHam, DiagonalMatrixElements) {
    const uint_t nwave = 3;
    InteractingBoseGasBosHam bos_ham(1, nwave, 1.0);
    Hamiltonian ham_src(&bos_ham);
    buffered::BosOnv mbf(ham_src.m_basis);
    mbf = {0, 0, 0, 0, 0, 3, 0};
    ham_t helem;
    helem = ham_src.get_element(mbf);
    ASSERT_NEAR_EQ(helem, 12);
}