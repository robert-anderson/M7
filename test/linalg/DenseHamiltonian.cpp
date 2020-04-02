//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include <src/linalg/EigenSolver.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/linalg/DenseHamiltonian.h"

#if 0
TEST(DenseHamiltonian, FciEnergyCheck) {
    DenseHamiltonian ham(AbInitioHamiltonian(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP"));
    auto solver = ham.diagonalize();
    // compare the ground and first excited states against BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals(0), -14.405976034322773));
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals(1), -14.288836984063053));
}
#endif