//
// Created by rja on 26/08/2021.
//

#include "src/core/excititer/ExcitIters.h"
#include "gtest/gtest.h"
#include "ExcitGenTester.h"
#include "src/core/excitgen/FrmBosKinetic.h"

TEST(FrmBosKinetic, HubbardKinetic){
    PRNG prng(14, 1000000);
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_KINETIC";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_NULL";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 3);
    excititers::FrmBosKinetic excit_iter(ham, exsig_utils::ex_1101);
    FrmBosKinetic excit_gen(ham, prng);
    excit_gen_tester::ExcitGenTester tester(excit_gen, excit_iter);
    buffered::FrmBosOnv src(ham.nsite());
    /*
     * arbitrary source ONV
     */
    src = {{0, 1, 4, 5}, {1, 0, 2}};
    tester.fill_results_table(src);
    std::cout << tester.m_results.to_string() << std::endl;
    tester.run(src, 100000);
    std::cout << tester.m_results.to_string() << std::endl;
}
