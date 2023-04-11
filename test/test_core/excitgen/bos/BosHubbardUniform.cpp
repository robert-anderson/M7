//
// Created by rja on 28/07/22.
//

#include "gtest/gtest.h"
#include "M7_lib/excitgen/bos/BosHubbardUniform.h"
#include "test_core/excitgen/ExcitGenTester.h"

TEST(BosHubbardUniform, Pbc1D) {
    PRNG prng(14, 1000000);
    HubbardBosHam bos_ham(4.0, lattice::make("ortho", {6}, {1}));
    Hamiltonian h(&bos_ham);
    exgen::BosHubbardUniform excit_gen(h.m_bos, prng);

    conn_foreach::bos::Hubbard conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::BosOnv src_mbf(h.m_basis);
    src_mbf = {1, 0, 2, 1, 3, 1};
    const auto correct = mpi::all_land(tester.run(src_mbf, 3000000) == "");
    ASSERT_TRUE(correct);
}