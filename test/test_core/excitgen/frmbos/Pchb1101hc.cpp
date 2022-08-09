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
    const uint_t nsite = 6;
    sys::frm::Basis frm_basis(lattice::make("ortho", {nsite}, {0}));
    sys::bos::Basis bos_basis(frm_basis.m_nsite, 4);
    GeneralLadderHam frmbos_ham({frm_basis, bos_basis}, {PROJECT_ROOT"/assets/SpinResolvedEbdump/EBDUMP", false});
    Hamiltonian h(&frmbos_ham);
    auto particles = h.default_particles(6);
    ASSERT_TRUE(frmbos_ham.m_basis.m_frm.m_abgrp_map);
    exgen::Pchb1101hc excit_gen(h.m_frmbos, prng);
    buffered::FrmBosOnv src_mbf(h.m_basis);
    mbf::set_aufbau_mbf(src_mbf.m_frm, particles.m_frm);
    ASSERT_EQ(src_mbf.m_frm.nsetbit(), particles.m_frm);
    typedef conn_foreach::frm::Ms2Conserve<1> frm_foreach_t;
    typedef conn_foreach::bos::Cre bos_foreach_t;
    typedef conn_foreach::frmbos::Product<frm_foreach_t, bos_foreach_t> foreach_t;

    foreach_t conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    ASSERT_EQ(tester.run(src_mbf, 5000000), "");
}
