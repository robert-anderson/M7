//
// Created by rja on 26/08/2021.
//

#include "src/core/excititer/ExcitIters.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "src/core/excitgen/HubbardSingles.h"
#include "src/core/field/Mbf.h"

TEST(Hubbard1dSingles, ObcFromNeel) {
    PRNG prng(14, 1000000);
    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", false);
    HubbardSingles excit_gen(ham, prng);
    ASSERT_FALSE(excit_gen.m_pbc);
    excititers::Frm excit_iter(ham, exsig_utils::ex_single);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.m_bd);
    Sector sector{ham.nelec(), true, 0, 0, 0};
    mbf::set_neel_mbf(src_mbf, sector);
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

TEST(Hubbard1dSingles, PbcFromNeel) {
    PRNG prng(14, 1000000);
    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_6site_pbc/FCIDUMP", false);
    HubbardSingles excit_gen(ham, prng);
    ASSERT_TRUE(excit_gen.m_pbc);
    excititers::Frm excit_iter(ham, exsig_utils::ex_single);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(ham.m_bd);
    Sector sector{ham.nelec(), true, 0, 0, 0};
    mbf::set_neel_mbf(src_mbf, sector);
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
