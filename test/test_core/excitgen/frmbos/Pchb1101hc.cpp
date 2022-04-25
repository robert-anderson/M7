//
// Created by rja on 05/04/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frmbos/Pchb1101hc.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "test_core/excitgen/ExcitGenTester.h"

TEST(Pchb1101hc, Test){
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_nelec = 6;
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/SpinResolvedEbdump/EBDUMP";
    opts.verify();
    Hamiltonian h(opts);
    ASSERT_TRUE(h.m_frmbos.is<GeneralLadderHam>());
    ASSERT_FALSE(h.m_hs.m_frm.m_abgrp_map.m_site_irreps.empty());
    Pchb1101hc excit_gen(h.m_frmbos, prng);
    buffered::FrmBosOnv src_mbf(h.m_hs);
    mbf::set_aufbau_mbf(src_mbf.m_frm);
    std::cout << src_mbf << std::endl;
    ASSERT_EQ(src_mbf.m_frm.nsetbit(), h.m_hs.m_frm.m_nelec);
    typedef conn_foreach::frm::Ms2Conserve<1> frm_foreach_t;
    typedef conn_foreach::bos::Cre bos_foreach_t;
    typedef conn_foreach::frmbos::Product<frm_foreach_t, bos_foreach_t> foreach_t;

    foreach_t conn_iter(h.m_hs.m_frm.m_sites, h.m_hs.m_bos);
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