//
// Created by rja on 26/08/2021.
//

#include "gtest/gtest.h"
#include "test_core/excitgen2/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen2/frm/HubbardUniform2.h"

TEST(HubbardUniform, ObcFromNeel1D) {
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4.0;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian h(opts);
    HubbardUniform2 excit_gen(*h.m_frm, prng);
    ASSERT_FALSE(excit_gen.h_cast()->m_bcs[0]);
    conn_foreach::frm::Hubbard conn_iter(excit_gen.h_cast()->m_lattice);

    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_bd);
    mbf::set_neel_mbf(src_mbf, h);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 3000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}

TEST(HubbardUniform, PbcFromNeel2D) {
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4.0;
    opts.m_fermion.m_hubbard.m_site_shape = {3, 3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {-1, 1};
    opts.verify();
    Hamiltonian h(opts);
    HubbardUniform2 excit_gen(*h.m_frm, prng);
    conn_foreach::frm::Hubbard conn_iter(excit_gen.h_cast()->m_lattice);

    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_bd);
    mbf::set_neel_mbf(src_mbf, h);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 3000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}