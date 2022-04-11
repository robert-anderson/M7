//
// Created by rja on 04/04/2022.
//

#include <M7_lib/excitgen2/frm/UniformSingles2.h>
#include "gtest/gtest.h"
#include "test_core/excitgen2/ExcitGenTester.h"

TEST(UniformSingles, FromExcited){
    PRNG prng(14, 1000000);
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/H2O_RHF/FCIDUMP";
    opts.verify();
    Hamiltonian h(opts);
    UniformSingles2 excit_gen(*h.m_frm, prng);
    conn_foreach::frm::Ms2Conserve<1> excit_iter(h.m_bd.m_frm.m_nsite);
    excit_gen_tester::ExcitGenTester tester(h, excit_gen, excit_iter);
    buffered::FrmOnv src_mbf(h.m_bd);
    src_mbf = {{0, 1, 2, 3, 4}, {0, 1, 2, 3, 5}};
    tester.fill_results_table(src_mbf);
    tester.run(src_mbf, 10000000);
    ASSERT_TRUE(tester.all_drawn_at_least_once());
}
