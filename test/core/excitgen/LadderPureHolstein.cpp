//
// Created by rja on 16/12/2021.
//

#include <src/core/excitgen/LadderPureHolstein.h>
#include <src/core/excititer/LadderPureHolstein.h>
#include "gtest/gtest.h"
#include "ExcitGenTester.h"


TEST(LadderPureHolstein, HolsteinAnn){
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {1};
    opts.m_ladder.m_holstein_coupling = 0.3;
    opts.m_ladder.m_nboson_max = 4;
    opts.m_boson.m_holstein_omega = 1.4;
    opts.verify();
    Hamiltonian ham(opts);
    excititers::LadderPureHolstein excit_iter(ham, exsig_utils::ex_0001);
    LadderHolsteinAnn excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 2, 6, 10, 11}, {1, 3, 2, 1, 0, 1}};
    tester.fill_results_table(src);
    const size_t ndraw = 100000;
    tester.run(src, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(LadderPureHolstein, HolsteinCre){
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {1};
    opts.m_ladder.m_holstein_coupling = 0.3;
    opts.m_ladder.m_nboson_max = 4;
    opts.m_boson.m_holstein_omega = 1.4;
    opts.verify();
    Hamiltonian ham(opts);
    excititers::LadderPureHolstein excit_iter(ham, exsig_utils::ex_0010);
    LadderHolsteinCre excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 2, 6, 10, 11}, {1, 3, 2, 1, 0, 1}};
    tester.fill_results_table(src);
    const size_t ndraw = 100000;
    tester.run(src, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}