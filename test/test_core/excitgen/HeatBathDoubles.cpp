//
// Created by rja on 09/05/2020.
//

#include <src/core/excititer/ExcitIters.h>
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "src/core/excitgen/HeatBathDoubles.h"
#include "src/core/field/Mbf.h"

TEST(HeatBathDoubles, SmallFromHFDeterminant){
    PRNG prng(14, 1000000);
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_LiH_STO-3G/FCIDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    HeatBathDoubles excit_gen(ham, prng);
    excititers::Frm excit_iter(ham, exsig_utils::ex_double);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.m_bd);
    mbf::set_aufbau_mbf(src_mbf, ham);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 10000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    std::cout << tester.m_results.to_string() << std::endl;
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2*ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2*ndraw));
}

TEST(HeatBathDoubles, LargeFromHFDeterminant){
    PRNG prng(14, 1000000);
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_N2_CCPVDZ/FCIDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    HeatBathDoubles excit_gen(ham, prng);
    excititers::Frm excit_iter(ham, exsig_utils::ex_double);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.m_bd);
    mbf::set_aufbau_mbf(src_mbf, ham);
    tester.fill_results_table(src_mbf);
    tester.run(src_mbf, 50000000);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
}