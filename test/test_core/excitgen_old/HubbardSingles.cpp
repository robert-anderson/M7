//
// Created by Robert J. Anderson on 26/08/2021.
//

#include "M7_lib/excititer/ExcitIters.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "M7_lib/excitgen/HubbardUniform.h"
#include "M7_lib/field/Mbf.h"

TEST(HubbardUniform, ObcFromNeel1D) {
    PRNG prng(14, 1000000);
    conf::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4.0;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian ham(opts);
    HubbardUniform excit_gen(ham, prng);
    ASSERT_FALSE(excit_gen.h_cast()->m_bcs[0]);
    excititers::Frm excit_iter(ham, exsig_utils::ex_single);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.m_bd);
    mbf::set_neel_mbf(src_mbf, ham);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 3000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(HubbardUniform, PbcFromNeel1D) {
    PRNG prng(14, 1000000);
    conf::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4.0;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {1};
    opts.verify();
    Hamiltonian ham(opts);
    HubbardUniform excit_gen(ham, prng);
    ASSERT_TRUE(excit_gen.h_cast()->m_bcs[0]);
    excititers::Frm excit_iter(ham, exsig_utils::ex_single);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.m_bd);
    mbf::set_neel_mbf(src_mbf, ham);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 3000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}
