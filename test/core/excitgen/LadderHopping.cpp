//
// Created by rja on 26/08/2021.
//

#include "src/core/excititer/ExcitIters.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "src/core/excitgen/LadderHoppingPc.h"

TEST(LadderHopping, HubbardUniform1101){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 2);
    excititers::LadderHopping excit_iter(ham, exsig_utils::ex_1101);
    LadderHoppingUniform excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 0, 2}};
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

TEST(LadderHopping, HubbardUniform1110){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 2);
    excititers::LadderHopping excit_iter(ham, exsig_utils::ex_1110);
    LadderHoppingUniform excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 0, 2}};
    tester.fill_results_table(src);
    const size_t ndraw = 200000;
    tester.run(src, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(LadderHopping, HubbardPc1101){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 2);
    excititers::LadderHopping excit_iter(ham, exsig_utils::ex_1101);
    LadderHoppingPc excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 0, 2}};
    tester.fill_results_table(src);
    const size_t ndraw = 500000;
    tester.run(src, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src, ndraw);
    std::cout << tester.m_results.to_string() << std::endl;
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(LadderHopping, HubbardPc1110){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 2);
    excititers::LadderHopping excit_iter(ham, exsig_utils::ex_1110);
    LadderHoppingPc excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 0, 2}};
    tester.fill_results_table(src);
    const size_t ndraw = 700000;
    tester.run(src, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src, ndraw);
    std::cout << tester.m_results.to_string() << std::endl;
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}
