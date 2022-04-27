//
// Created by Robert J. Anderson on 28/08/2021.
//

#include <M7_lib/excitgen/LadderPureHolstein.h>
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/excititer/LadderPure.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"

TEST(LadderPure, Uniform0001){
    PRNG prng(14, 1000000);
    conf::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_hamiltonian.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HH_V1.4_WITH_UNC";
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    std::vector<defs::ham_t> uncs_chk = {0.4, 1.3, 0.7};
    auto& h_cast = dynamic_cast<const GeneralLadderHam&>(*ham.m_ladder);
    ASSERT_EQ(h_cast.m_v_unc, uncs_chk);
    excititers::LadderPure excit_iter(ham, exsig_utils::ex_0001);
    LadderPureUniform excit_gen(ham, prng, {exsig_utils::ex_0001});
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
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

TEST(LadderPure, Uniform0010){
    PRNG prng(14, 1000000);
    conf::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_hamiltonian.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HH_V1.4_WITH_UNC";
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    std::vector<defs::ham_t> uncs_chk = {0.4, 1.3, 0.7};
    auto& h_cast = dynamic_cast<const GeneralLadderHam&>(*ham.m_ladder);
    ASSERT_EQ(h_cast.m_v_unc, uncs_chk);
    excititers::LadderPure excit_iter(ham, exsig_utils::ex_0010);
    LadderPureUniform excit_gen(ham, prng, {exsig_utils::ex_0010});
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.m_bd);
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