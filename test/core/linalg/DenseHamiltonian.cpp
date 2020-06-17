//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include <src/core/linalg/EigenSolver.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/linalg/DenseHamiltonian.h"

TEST(DenseHamiltonian, FciEnergyCheck4c) {
    if (consts::is_complex<defs::ham_comp_t>()) {
        DenseHamiltonian ham(AbInitioHamiltonian(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP"));
        auto solver = ham.diagonalize();
        // compare the ground and first excited states to BAGEL's values
        ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -14.40597603432, 1e-10));
        ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[1], -14.28883698406, 1e-10));
    }
}

TEST(DenseHamiltonian, FciEnergyCheckRhf) {
    AbInitioHamiltonian ham_src(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP");
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -108.81138657563143, 1e-10));
}