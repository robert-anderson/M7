//
// Created by rja on 05/04/2022.
//

#include "gtest/gtest.h"
#include "test_core/excitgen2/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen2/frmbos/Pchb1101hc.h"

TEST(Pchb1101hc, Test){
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/SpinResolvedEbdump/EBDUMP";
    opts.verify();
    Hamiltonian h(opts);
    Pchb1101hc excit_gen(*h.m_frmbos, prng);
    buffered::FrmOnv src_mbf(h.m_bd);
    mbf::set_aufbau_mbf(src_mbf, h);
    conn_foreach::frm::Ms2Conserve<2> conn_iter(src_mbf.m_bd);
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);

    tester.fill_results_table(src_mbf);
    const size_t ndraw = 10000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}