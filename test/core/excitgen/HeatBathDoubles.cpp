//
// Created by rja on 09/05/2020.
//

#include <src/core/excititer/ExcitIters.h>
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "src/core/excitgen/HeatBathDoubles.h"

TEST(HeatBathDoubles, SmallFromHFDeterminant){
    PRNG prng(14, 1000000);
    Hamiltonian ham(defs::assets_root + "/RHF_LiH_STO-3G/FCIDUMP", false);
    HeatBathDoubles excit_gen(ham, prng);
    excititers::FrmConserve excit_iter(ham, exsig_utils::ex_double);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.nsite());
    ham.set_hf_mbf(src_mbf, 0);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 10000000;
    tester.run<field::FrmOnv>(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run<field::FrmOnv>(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2*ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2*ndraw));
}

TEST(HeatBathDoubles, LargeFromHFDeterminant){
    PRNG prng(14, 1000000);
    Hamiltonian ham(defs::assets_root + "/RHF_N2_CCPVDZ/FCIDUMP", false);
    HeatBathDoubles excit_gen(ham, prng);
    excititers::FrmConserve excit_iter(ham, exsig_utils::ex_double);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.nsite());
    ham.set_hf_mbf(src_mbf, 0);
    tester.fill_results_table(src_mbf);
    tester.run<field::FrmOnv>(src_mbf, 10000000);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
}