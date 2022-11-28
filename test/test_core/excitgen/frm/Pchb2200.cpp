//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "gtest/gtest.h"
#include "test_core/excitgen/ExcitGenTester.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/excitgen/frm/Pchb2200.h"

TEST(Pchb2200, FromHFDeterminant) {
    PRNG prng(14, 1000000);
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_LiH_STO-3G/FCIDUMP"});
    Hamiltonian h(&frm_ham);
    exgen::Pchb2200 excit_gen(frm_ham, prng);
    conn_foreach::frm::Ms2Conserve<2> conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(frm_ham.m_basis);
    mbf::set_aufbau_mbf(src_mbf, h.default_particles().m_frm);
    ASSERT_EQ(tester.run(src_mbf, 10000000), "");
}

TEST(Pchb2200, FromExcited){
    PRNG prng(14, 1000000);
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    Hamiltonian h(&frm_ham);
    exgen::Pchb2200 excit_gen(frm_ham, prng);
    conn_foreach::frm::Ms2Conserve<2> conn_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, conn_iter);
    buffered::FrmOnv src_mbf(frm_ham.m_basis);
    src_mbf = {{0, 1, 3}, {1, 2, 4}};
    ASSERT_EQ(tester.run(src_mbf, 50000000), "");
}
