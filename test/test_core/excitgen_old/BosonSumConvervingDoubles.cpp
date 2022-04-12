//
// Created by rja on 26/11/2021.
//

#include <M7_lib/excititer/ExcitIters.h>
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "M7_lib/excitgen/BosonSumConservingDoubles.h"

TEST(BosonSumConservingDoubles, LandauLevels) {
    PRNG prng(14, 1000000);
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/LandauLevels_6_8_14/BOSDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    BosonSumConservingDoubles excit_gen(ham, prng);
    excititers::Bos excit_iter(ham, exsig_utils::ex_0022);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::BosOnv src_mbf(ham.m_bd);
    src_mbf = {2, 0, 1, 2, 0, 0, 1, 0};
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 10000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    std::cout << tester.m_results.to_string() << std::endl;
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LT(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}
