//
// Created by Robert J. Anderson on 26/08/2021.
//

#include "M7_lib/excititer/LadderHopping.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "M7_lib/excitgen/LadderHoppingPc.h"

TEST(LadderHopping, HubbardUniform1101){
    PRNG prng(14, 1000000);
    conf::Document opts;
    opts.m_hamiltonian.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_hamiltonian.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_hamiltonian.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_hamiltonian.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
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
    conf::Document opts;
    opts.m_hamiltonian.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_hamiltonian.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_hamiltonian.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_hamiltonian.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
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
    conf::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_hamiltonian.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
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
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(LadderHopping, HubbardPc1110){
    PRNG prng(14, 1000000);
    conf::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_hamiltonian.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HOPPING";
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
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
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}
