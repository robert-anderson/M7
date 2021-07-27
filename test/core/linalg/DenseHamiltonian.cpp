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

TEST(DenseHamiltonian, N2Rhf) {
    Hamiltonian ham_src(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -108.916561245585, 1e-8));
}

TEST(DenseHamiltonian, HF) {
    Hamiltonian ham_src(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -99.9421389039332, 1e-8));
}

TEST(DenseHamiltonian, PyscfX2cCheck) {
    Hamiltonian ham_src(defs::assets_root + "/H2O_X2C/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -76.08150945314577, 1e-10));
}

TEST(DenseHamiltonian, Hubbard3Site) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_3site/FCIDUMP", 1);
    ASSERT_EQ(h.nelec(), 4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], 2.0);
}

TEST(DenseHamiltonian, Hubbard4Site) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 1);
    ASSERT_EQ(h.nelec(), 4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, Hubbard6Site) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 1);
    ASSERT_EQ(h.nelec(), 6);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -3.0925653194551845);
}

TEST(DenseHamiltonian, BosonCouplingNoBosonLimit) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 0);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, BosonCouplingNoField) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_3site/FCIDUMP", 0, 2);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], 2.0);
}

TEST(DenseHamiltonian, BosonCouplingNoFrequencyMaxOcc1) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 1, 0.0, 1.4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -7.553145308678367);
}

TEST(DenseHamiltonian, BosonCouplingNoFrequencyMaxOcc2) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_3site/FCIDUMP", 0, 2, 0.0, 1.4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -7.699484522379835);
}

TEST(DenseHamiltonian, BosonCouplingNoFrequencyMaxOcc3) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_3site/FCIDUMP", 0, 3, 0.0, 1.4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -11.07271962268484);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc1) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 1, 0.3, 1.4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -6.9875779675355165);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc2) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_3site/FCIDUMP", 0, 2, 0.3, 1.4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -6.692966463435127);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc3) {
    Hamiltonian h(defs::assets_root + "/Hubbard_U4_3site/FCIDUMP", 0, 3, 0.3, 1.4);
    DenseHamiltonian dh(h);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -9.423844225360671);
}