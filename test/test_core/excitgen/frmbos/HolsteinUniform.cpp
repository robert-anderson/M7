//
// Created by rja on 21/07/22.
//

#include "gtest/gtest.h"
#include "M7_lib/excitgen/frmbos/HolsteinUniform.h"
#include "M7_lib/hamiltonian/Hamiltonian.h"
#include "test_core/excitgen/ExcitGenTester.h"

TEST(HolsteinUniform, Creation) {
    PRNG prng(14, 1000000);
    uint_t nsite = 8;
    HolsteinLadderHam frmbos_ham({nsite}, 1.0, 4);
    Hamiltonian h(&frmbos_ham);
    exgen::HolsteinUniform0010 excit_gen(h.m_frmbos, prng);
    buffered::FrmBosOnv src_mbf(h.m_basis);
    src_mbf.m_frm = {{0, 1, 4, 7}, {2, 4, 5, 7}};
    src_mbf.m_bos = {0, 2, 3, 1, 0, 4, 1, 2};
    conn_foreach::bos::Cre conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    ASSERT_EQ(tester.run(src_mbf, 5000000), "");
}

TEST(HolsteinUniform, Annihilation) {
    PRNG prng(14, 1000000);
    uint_t nsite = 8;
    HolsteinLadderHam frmbos_ham({nsite}, 1.0, 4);
    Hamiltonian h(&frmbos_ham);
    exgen::HolsteinUniform0001 excit_gen(h.m_frmbos, prng);
    buffered::FrmBosOnv src_mbf(h.m_basis);
    src_mbf.m_frm = {{0, 1, 4, 7}, {2, 4, 5, 7}};
    src_mbf.m_bos = {0, 2, 3, 1, 0, 4, 1, 2};
    conn_foreach::bos::Ann conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    ASSERT_EQ(tester.run(src_mbf, 5000000), "");
}