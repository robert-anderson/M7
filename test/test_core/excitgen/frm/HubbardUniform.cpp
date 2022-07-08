//
// Created by Robert J. Anderson on 26/08/2021.
//

#include "gtest/gtest.h"
#include "test_core/excitgen/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frm/HubbardUniform.h"

TEST(HubbardUniform, ObcFromNeel1D) {
    PRNG prng(14, 1000000);
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {6}, {0}));
    Hamiltonian h(&frm_ham);
    exgen::HubbardUniform excit_gen(h.m_frm, prng);
    ASSERT_FALSE(frm_ham.m_basis.m_lattice->as<lattice::Ortho>()->m_topo.m_bcs[0]);

    conn_foreach::frm::Hubbard conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_basis);
    mbf::set_neel_mbf(src_mbf, h.default_particles().m_frm);
    ASSERT_EQ(tester.run(src_mbf, 3000000), "");
}

TEST(HubbardUniform, PbcFromNeel2D) {
    PRNG prng(14, 1000000);
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3, 3}, {-1, 1}));
    Hamiltonian h(&frm_ham);
    exgen::HubbardUniform excit_gen(h.m_frm, prng);
    conn_foreach::frm::Hubbard conn_iter;

    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_basis);
    mbf::set_neel_mbf(src_mbf, h.default_particles().m_frm);
    ASSERT_EQ(tester.run(src_mbf, 3000000), "");
}

TEST(HubbardPreferDoubleOcc, Pbc2D) {
    PRNG prng(14, 1000000);
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3, 3}, {-1, 1}));
    Hamiltonian h(&frm_ham);
    exgen::HubbardPreferDoubleOcc excit_gen(h.m_frm, prng, 1.0);
    conn_foreach::frm::Hubbard conn_iter;

    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(h.m_basis);
    const uintv_t alpha_sites = {0, 1, 0, 0, 0, 1, 1, 1, 0};
    const uintv_t beta_sites =  {1, 1, 0, 1, 0, 0, 0, 1, 1};
    src_mbf.set(alpha_sites, beta_sites);
    mbf::set_neel_mbf(src_mbf, h.default_particles().m_frm);
    ASSERT_EQ(tester.run(src_mbf, 3000000), "");
}