//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include <src/core/linalg/EigenSolver.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/linalg/DenseHamiltonian.h"

TEST(DenseHamiltonian, FciEnergyCheck) {
    DenseHamiltonian ham(AbInitioHamiltonian(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP"));
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    std::cout << std::endl;
    std::cout << solver.m_evals(0) <<std::endl;
    std::cout << solver.m_evals(1) <<std::endl;
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals(0), -14.405976034322773));
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals(1), -14.288836984063053));
}