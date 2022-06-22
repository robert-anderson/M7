//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "gtest/gtest.h"
#include "test_core/excitgen/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frm/Pchb2200.h"

TEST(Pchb2200, FromHFDeterminant) {
    PRNG prng(14, 1000000);
    GeneralFrmHam frm_ham({defs::c_assets_root + "/RHF_LiH_STO-3G/FCIDUMP"}, true);
    Hamiltonian h(&frm_ham);
    Pchb2200 excit_gen(frm_ham, prng);
    conn_foreach::frm::Ms2Conserve<2> conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(frm_ham.m_basis);
    mbf::set_aufbau_mbf(src_mbf, h.default_particles().m_frm);

    tester.fill_results_table(src_mbf);
    const size_t ndraw = 10000000;
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(Pchb2200, FromExcited){
    PRNG prng(14, 1000000);
    GeneralFrmHam frm_ham({defs::c_assets_root + "/RHF_N2_6o6e/FCIDUMP"}, true);
    Hamiltonian h(&frm_ham);
    Pchb2200 excit_gen(frm_ham, prng);
    conn_foreach::frm::Ms2Conserve<2> conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(frm_ham.m_basis);
    src_mbf = {{0, 1, 3}, {1, 2, 4}};
    tester.fill_results_table(src_mbf);
    ASSERT_EQ(tester.run(src_mbf, 50000000).m_error_message, "");
    ASSERT_TRUE(tester.all_drawn_at_least_once());
}
