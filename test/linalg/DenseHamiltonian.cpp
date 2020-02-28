//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include "src/hamiltonian/AbInitioHamiltonian.h"
#include "src/linalg/DenseHamiltonian.h"


TEST(DenseHamiltonian, FciEnergyCheck) {
    DenseHamiltonian ham(AbInitioHamiltonian(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP"));
    auto es = ham.diagonalise();
    // compare the ground and first excited states against BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(es.get()->eigenvalues()[0], -14.405976034322808));
    ASSERT_TRUE(consts::floats_nearly_equal(es.get()->eigenvalues()[1], -14.288836984063064));
}