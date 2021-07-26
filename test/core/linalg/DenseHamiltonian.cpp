//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include <src/core/linalg/EigenSolver.h>
#include "src/core/linalg/DenseHamiltonian.h"

#ifdef ENABLE_COMPLEX
TEST(DenseHamiltonian, FciEnergyCheck4c) {
    DenseHamiltonian ham(FermionHamiltonian(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false));
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -14.40597603432, 1e-10));
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[1], -14.28883698406, 1e-10));
}
#endif

TEST(DenseHamiltonian, FciEnergyCheckRhf) {
    FermionHamiltonian ham_src(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -108.916561245585, 1e-8));
}

TEST(DenseHamiltonian, FciEnergyCheckHf) {
    FermionHamiltonian ham_src(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -99.9421389039332, 1e-8));
}

TEST(DenseHamiltonian, PyscfX2cCheck) {
    FermionHamiltonian ham_src(defs::assets_root + "/H2O_X2C/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -76.08150945314577, 1e-10));
}

TEST(DenseHamiltonian, HubbardCheck) {
    FermionHamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 1);
    ASSERT_EQ(h.nelec(), 4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, HubbardCheck6Site) {
    FermionHamiltonian h(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 1);
    ASSERT_EQ(h.nelec(), 6);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    std::cout << solver.m_evals[0] << std::endl;
    ASSERT_FLOAT_EQ(solver.m_evals[0], -3.0925653194551845);
}

#ifdef ENABLE_BOSONS
TEST(DenseHamiltonian, BosonCouplingNoBosonLimit) {
    FermiBosHamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 0, 0, 0);
    DenseHamiltonian dh(h, 0);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, BosonCouplingNoField) {
    FermiBosHamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 2, 0, 0);
    DenseHamiltonian dh(h, 0);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc1) {
    FermiBosHamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 1, 1.4, 0.3);
    DenseHamiltonian dh(h, 0);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -6.9875779675355165);
}

TEST(DenseHamiltonian, BosonCouplingNoFrequencyMaxOcc2) {
    FermiBosHamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 2, 1.4, 0.0);
    DenseHamiltonian dh(h, 0);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -11.652629830979253);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc2) {
    FermiBosHamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 2, 1.4, 0.3);
    DenseHamiltonian dh(h, 0);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -10.328242246088791);
}
#endif