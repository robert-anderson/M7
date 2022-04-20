//
// Created by rja on 05/04/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frm/HeisenbergExchange.h"
#include "test_core/excitgen/ExcitGenTester.h"

TEST(HeisenbergExchange, Pbc2D) {
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_heisenberg.m_site_shape = {3, 3};
    opts.m_fermion.m_heisenberg.m_boundary_conds = {1, 1};
    opts.m_fermion.m_ms2_restrict = 1;
    opts.verify();
    Hamiltonian h(opts);
    ASSERT_TRUE(dynamic_cast<const HeisenbergFrmHam*>(&h.m_frm));
    HeisenbergExchange excit_gen(h.m_frm, prng);
    conn_foreach::frm::Heisenberg conn_iter(excit_gen.h_cast()->m_lattice);

    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_bd);
    mbf::set_neel_mbf(src_mbf, h);
    tester.fill_results_table(src_mbf);
    const size_t ndraw = 1000000;
    tester.run(src_mbf, ndraw);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    tester.run(src_mbf, ndraw);
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}
