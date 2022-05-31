//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "gtest/gtest.h"

#include "M7_lib/excitgen/frm/UniformSingles.h"
#include "test_core/excitgen/ExcitGenTester.h"

TEST(UniformSingles, FromExcited){
    PRNG prng(14, 1000000);
    GeneralFrmHam frm_ham({defs::assets_root + "/H2O_RHF/FCIDUMP"}, true);
    Hamiltonian h(&frm_ham);
    UniformSingles excit_gen(frm_ham, prng);
    conn_foreach::frm::Ms2Conserve<1> excit_iter;
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(h.m_basis);
    src_mbf = {{0, 1, 2, 3, 4}, {0, 1, 2, 3, 5}};
    tester.fill_results_table(src_mbf);
    tester.run(src_mbf, 10000000);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
}
