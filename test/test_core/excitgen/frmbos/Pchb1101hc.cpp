//
// Created by Robert J. Anderson on 05/04/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frmbos/Pchb1101hc.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "test_core/excitgen/ExcitGenTester.h"

TEST(Pchb1101hc, Test){
    PRNG prng(14, 1000000);
    GeneralLadderHam frmbos_ham({defs::assets_root + "/SpinResolvedEbdump/EBDUMP"}, false);
    Hamiltonian h(&frmbos_ham);
    auto particles = h.default_particles(6);
    ASSERT_TRUE(frmbos_ham.m_basis.m_frm.m_abgrp_map);
    Pchb1101hc excit_gen(h.m_frmbos, prng);
    buffered::FrmBosOnv src_mbf(h.m_basis);
    mbf::set_aufbau_mbf(src_mbf.m_frm, particles.m_frm);
    ASSERT_EQ(src_mbf.m_frm.nsetbit(), particles.m_frm);
    typedef conn_foreach::frm::Ms2Conserve<1> frm_foreach_t;
    typedef conn_foreach::bos::Cre bos_foreach_t;
    typedef conn_foreach::frmbos::Product<frm_foreach_t, bos_foreach_t> foreach_t;

    foreach_t conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);

    tester.fill_results_table(src_mbf);
    const size_t ndraw = 5000000;
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    ASSERT_TRUE(tester.all_drawn_at_least_once());
    auto av_err1 = tester.mean_abs_error(ndraw);
    ASSERT_EQ(tester.run(src_mbf, ndraw).m_error_message, "");
    auto av_err2 = tester.mean_abs_error(2 * ndraw);
    ASSERT_LE(av_err2, av_err1);
    ASSERT_TRUE(tester.all_correct_weights(2 * ndraw));
}
