//
// Created by rja on 28/08/2021.
//

#include <src/core/excitgen/LadderPureHolsteinZpm.h>
#include "src/core/excititer/ExcitIters.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "src/core/excitgen/LadderPureUniform.h"

TEST(LadderPure, HubbardUniform0001){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HH_V1.4_WITH_UNC";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 3);
    std::vector<defs::ham_t> uncs_chk = {0.4, 1.3, 0.7};
    ASSERT_EQ(ham.m_ladder.m_v_unc, uncs_chk);
    excititers::LadderPure excit_iter(ham, exsig_utils::ex_0001);
    LadderPureUniform excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.nsite());
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 1, 2}};
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

TEST(LadderPure, HubbardUniform0010){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HH_V1.4_WITH_UNC";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 3);
    excititers::LadderPure excit_iter(ham, exsig_utils::ex_0010);
    LadderPureUniform excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.nsite());
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 1, 2}};
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

TEST(LadderPure, HolsteinZpm0001){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/HH_ZPM/FCIDUMP_ZPM";
    auto fname_eb = defs::assets_root + "/HH_ZPM/EBDUMP_ZPM";
    auto fname_bos = defs::assets_root + "/HH_ZPM/BOSDUMP";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 3);
    std::vector<defs::ham_t> uncs_chk = {-0.1, -0.1, -0.1, -0.1};
    ASSERT_EQ(ham.m_ladder.m_v_unc, uncs_chk);
    excititers::LadderPure excit_iter(ham, exsig_utils::ex_0001);
    LadderPureHolsteinZpm excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.nsite());
    /*
     * arbitrary source ONV
     */
    src = {{0, 3, 4, 5, 7}, {2, 1, 1, 1}};
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

TEST(LadderPure, HolsteinZpm0010){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/HH_ZPM/FCIDUMP_ZPM";
    auto fname_eb = defs::assets_root + "/HH_ZPM/EBDUMP_ZPM";
    auto fname_bos = defs::assets_root + "/HH_ZPM/BOSDUMP";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 3);
    std::vector<defs::ham_t> uncs_chk = {-0.1, -0.1, -0.1, -0.1};
    ASSERT_EQ(ham.m_ladder.m_v_unc, uncs_chk);
    excititers::LadderPure excit_iter(ham, exsig_utils::ex_0010);
    LadderPureHolsteinZpm excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.nsite());
    /*
     * arbitrary source ONV
     */
    src = {{0, 3, 4, 5, 7}, {2, 1, 1, 1}};
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