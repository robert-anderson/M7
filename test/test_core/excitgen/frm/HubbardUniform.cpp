//
// Created by Robert J. Anderson on 26/08/2021.
//

#include "gtest/gtest.h"
#include "test_core/excitgen/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frm/HubbardUniform.h"

TEST(HubbardUniform, ObcFromNeel1D) {
    PRNG prng(14, 1000000);
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {6}, {0}));
    Hamiltonian h(&frm_ham);
    HubbardUniform excit_gen(h.m_frm, prng);
    ASSERT_FALSE(frm_ham.m_basis.m_lattice->as<lattice::Ortho>()->m_topo.m_bcs[0]);

    conn_foreach::frm::Hubbard conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_basis);
    mbf::set_neel_mbf(src_mbf, h.default_particles().m_frm);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 3000000;
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(HubbardUniform, PbcFromNeel2D) {

    PRNG prng(14, 1000000);
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3, 3}, {-1, 1}));
    Hamiltonian h(&frm_ham);
    HubbardUniform excit_gen(h.m_frm, prng);
    conn_foreach::frm::Hubbard conn_iter;

    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_basis);
    mbf::set_neel_mbf(src_mbf, h.default_particles().m_frm);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 3000000;
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}