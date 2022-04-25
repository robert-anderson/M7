//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include "M7_lib/linalg/DenseHamiltonian.h"

/*
 * exact diagonalization in the entire Hilbert space for integration testing of matrix elements for very small systems
 */

#ifdef ENABLE_COMPLEX
TEST(DenseHamiltonian, FciEnergyCheck4c) {
    DenseHamiltonian ham(Hamiltonian(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false));
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -14.40597603432, 1e-10));
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[1], -14.28883698406, 1e-10));
}
#endif
TEST(DenseHamiltonian, N2Rhf) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_N2_6o6e/FCIDUMP";
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<defs::ham_t> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -108.916561245585, 1e-8));
}

TEST(DenseHamiltonian, HeisenbergFrmHam) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_heisenberg.m_boundary_conds = {1};
    opts.verify();
    /*
     * https://doi.org/10.1016/0378-4363(78)90115-8
     */
    defs::inds nsites = {4, 6, 8, 10, 12};
    std::vector<defs::ham_comp_t> energies = {-2.0, -2.8027756375, -3.6510934085, -4.515446354, -5.387390917};
    for (size_t i=0ul; i<nsites.size(); ++i){
        auto nsite = nsites[i];
        auto energy = energies[i];
        opts.m_fermion.m_heisenberg.m_site_shape = {nsite};
        Hamiltonian ham_src(opts);
        std::vector<double> evals;
        DenseHamiltonian ham(ham_src);
        dense::diag(ham, evals);
        ASSERT_TRUE(consts::nearly_equal(evals[0], energy, 1e-8));
    }
}

TEST(DenseHamiltonian, HF) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/HF_RDMs/FCIDUMP";
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -99.9421389039332, 1e-8));
}

TEST(DenseHamiltonian, PyscfX2cCheck) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/H2O_X2C/FCIDUMP";
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::nearly_equal(evals[0], -76.08150945314577, 1e-10));
}

TEST(DenseHamiltonian, Hubbard3Site) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_fermion.m_nelec = -1;
    opts.verify();
    Hamiltonian ham_src(opts);
    ASSERT_EQ(ham_src.m_hs.m_frm.m_nelec, 4);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], 2.0, 1e-10));
}

TEST(DenseHamiltonian, Hubbard4Site) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {4};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian ham_src(opts);
    ASSERT_EQ(ham_src.m_hs.m_frm.m_nelec, 4);
    ASSERT_EQ(ham_src.m_hs.m_frm.m_sites, 4);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -1.9531453086749293, 1e-10));
}

TEST(DenseHamiltonian, Hubbard6Site) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian ham_src(opts);
    ASSERT_EQ(ham_src.m_hs.m_frm.m_nelec, 6);
    ASSERT_EQ(ham_src.m_hs.m_frm.m_sites, 6);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -3.0925653194551845, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinNoCoupling) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_nelec = -1;
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 0.0;
    opts.m_boson.m_bos_occ_cutoff = 0;
    opts.m_boson.m_num_op_weight = 0.0;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], 2.0, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinNoFrequencyMaxOcc2) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_nelec = -1;
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_boson.m_bos_occ_cutoff = 2;
    opts.m_boson.m_num_op_weight = 0.0;
    opts.verify();
    Hamiltonian ham_src(opts);
    const auto nmode = ham_src.m_hs.m_bos.m_nmode;
    const auto& sites = ham_src.m_hs.m_frm.m_sites;
    for (size_t n = 0ul; n < nmode; ++n) {
        for (size_t p = 0ul; p < sites.m_nspinorb; ++p) {
            const auto psite = sites.isite(p);
            for (size_t q = 0ul; q < sites.m_nspinorb; ++q) {
                const auto qsite = sites.isite(q);
                if (n == psite && psite == qsite) {
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_frmbos.get_coeff_1101(n, p, q)), 1.4);
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_frmbos.get_coeff_1110(n, p, q)), 1.4);
                } else {
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_frmbos.get_coeff_1101(n, p, q)), 0.0);
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_frmbos.get_coeff_1110(n, p, q)), 0.0);
                }
            }
        }
    }
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -7.699484522379835, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinNoFrequencyMaxOcc3) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_nelec = -1;
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_boson.m_bos_occ_cutoff = 3;
    opts.m_boson.m_num_op_weight = 0.0;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -11.07271962268484, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinMaxOcc2) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_nelec = -1;
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_boson.m_bos_occ_cutoff = 2;
    opts.m_boson.m_num_op_weight = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -6.692966463435127, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinMaxOcc1) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_nelec = -1;
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_boson.m_bos_occ_cutoff = 1;
    opts.m_boson.m_num_op_weight = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -3.1699561178752873, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinMaxOcc3) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_nelec = -1;
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_boson.m_bos_occ_cutoff = 3;
    opts.m_boson.m_num_op_weight = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -9.423844225360671, 1e-10));
}

TEST(DenseHamiltonian, BosonCouplingGeneralMaxOcc1) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL";
    opts.m_boson.m_bos_occ_cutoff = 1;
    opts.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL";
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], 0.5090148148366922, 1e-10));
}

TEST(DenseHamiltonian, BosonCouplingGeneralMaxOcc2) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL";
    opts.m_boson.m_bos_occ_cutoff = 2;
    opts.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL";
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -0.38125085248276913, 1e-10));
}

TEST(DenseHamiltonian, BosonCouplingGeneralMaxOcc3) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL";
    opts.m_boson.m_bos_occ_cutoff = 3;
    opts.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL";
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -0.9998830020871416, 1e-10));
}