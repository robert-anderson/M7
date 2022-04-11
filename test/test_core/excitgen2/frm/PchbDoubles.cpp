//
// Created by rja on 04/04/2022.
//

#include "gtest/gtest.h"
#include "test_core/excitgen2/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen2/frm/Pchb2200.h"

TEST(HeatBathDoubles, FromHFDeterminant) {
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_LiH_STO-3G/FCIDUMP";
    opts.verify();
    Hamiltonian h(opts);
    Pchb2200 excit_gen(*h.m_frm, prng);
    conn_foreach::frm::Ms2Conserve<2> conn_iter(10);
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_bd);
    mbf::set_aufbau_mbf(src_mbf, h);

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

TEST(HeatBathDoubles, FromExcited){
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_N2_6o6e/FCIDUMP";
    opts.verify();
    Hamiltonian h(opts);
    Pchb2200 excit_gen(*h.m_frm, prng);
    conn_foreach::frm::Ms2Conserve<2> excit_iter(h.m_bd.m_frm.m_nsite);
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(h.m_bd);
    src_mbf = {{0, 1, 3}, {1, 2, 4}};
    tester.fill_results_table(src_mbf);
    tester.run(src_mbf, 50000000);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
}
