//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include <src/core/linalg/EigenSolver.h>
#include "src/core/linalg/DenseHamiltonian.h"

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
    auto solver = ham.diagonalize();
    ASSERT_TRUE(consts::nearly_equal(solver.m_evals[0], -108.916561245585, 1e-8));
}

TEST(DenseHamiltonian, HeisenbergFrmHam) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_heisenberg.m_site_shape = {4};
    opts.m_fermion.m_heisenberg.m_boundary_conds = {1};
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    /*
     * https://doi.org/10.1016/0378-4363(78)90115-8
     *  nsite  exact eigenvalue
     *   4     -4.0
     *   6     -5.605551275463988
     *   8     -7.3021868178743485
     *   10    -9.030892708984082
     *   12    -10.774781834890415
     */
    std::cout << solver.m_evals[0] << std::endl;
    ASSERT_TRUE(consts::nearly_equal(solver.m_evals[0], -5.605551275463988, 1e-8));
}

#if 0
TEST(DenseHamiltonian, HF) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/HF_RDMs/FCIDUMP";
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -99.9421389039332, 1e-8));
}

TEST(DenseHamiltonian, PyscfX2cCheck) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/H2O_X2C/FCIDUMP";
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -76.08150945314577, 1e-10));
}

TEST(DenseHamiltonian, Hubbard3Site) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian ham_src(opts);
    ASSERT_EQ(ham_src.nelec(), 4);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], 2.0);
}

TEST(DenseHamiltonian, Hubbard4Site) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {4};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian ham_src(opts);
    ASSERT_EQ(ham_src.nelec(), 4);
    ASSERT_EQ(ham_src.m_bd.m_nsite, 4);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, Hubbard6Site) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {6};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.verify();
    Hamiltonian ham_src(opts);
    ASSERT_EQ(ham_src.nelec(), 6);
    ASSERT_EQ(ham_src.m_bd.m_nsite, 6);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -3.0925653194551845);
}

TEST(DenseHamiltonian, BosonCouplingNoField) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 0.0;
    opts.m_boson.m_holstein_omega = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], 2.0);
}

TEST(DenseHamiltonian, BosonCouplingNoFrequencyMaxOcc2) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_ladder.m_nboson_max = 2;
    opts.m_boson.m_holstein_omega = 0.0;
    opts.verify();
    Hamiltonian ham_src(opts);
    for (size_t n = 0ul; n < ham_src.m_bd.m_nmode; ++n) {
        for (size_t p = 0ul; p < ham_src.m_bd.m_nspinorb; ++p) {
            const auto psite = FrmOnvField::isite(p, ham_src.m_bd.m_nsite);
            for (size_t q = 0ul; q < ham_src.m_bd.m_nspinorb; ++q) {
                const auto qsite = FrmOnvField::isite(q, ham_src.m_bd.m_nsite);
                if (n == psite && psite == qsite) {
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_ladder->get_coeff_1101(n, p, q)), 1.4);
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_ladder->get_coeff_1110(n, p, q)), 1.4);
                } else {
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_ladder->get_coeff_1101(n, p, q)), 0.0);
                    ASSERT_FLOAT_EQ(consts::real(ham_src.m_ladder->get_coeff_1110(n, p, q)), 0.0);
                }
            }
        }
    }
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -7.699484522379835);
}

TEST(DenseHamiltonian, BosonCouplingNoFrequencyMaxOcc3) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_ladder.m_nboson_max = 3;
    opts.m_boson.m_holstein_omega = 0.0;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -11.07271962268484);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc2) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_ladder.m_nboson_max = 2;
    opts.m_boson.m_holstein_omega = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -6.692966463435127);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc1) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_ladder.m_nboson_max = 1;
    opts.m_boson.m_holstein_omega = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -3.1699561178752873);
}

TEST(DenseHamiltonian, BosonCouplingMaxOcc3) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_hubbard.m_repulsion = 4;
    opts.m_fermion.m_hubbard.m_site_shape = {3};
    opts.m_fermion.m_hubbard.m_boundary_conds = {0};
    opts.m_ladder.m_holstein_coupling = 1.4;
    opts.m_ladder.m_nboson_max = 3;
    opts.m_boson.m_holstein_omega = 0.3;
    opts.verify();
    Hamiltonian ham_src(opts);
    DenseHamiltonian dh(ham_src);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -9.423844225360671);
}

TEST(DenseHamiltonian, BosonCouplingGeneralMaxOcc1) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL";
    opts.m_ladder.m_nboson_max = 1;
    opts.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL";
    Hamiltonian src_ham(opts);
    DenseHamiltonian dh(src_ham);
    const auto frm_dim = ci_utils::fermion_dim(src_ham.m_bd.m_nsite, src_ham.nelec());
    const auto bos_dim = ci_utils::boson_dim(src_ham.m_bd.m_nsite, src_ham.m_nboson_max, false);
    ASSERT_EQ(dh.m_ncol, frm_dim * bos_dim);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], 0.5090148148366922);
}

TEST(DenseHamiltonian, BosonCouplingGeneralMaxOcc2) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL";
    opts.m_ladder.m_nboson_max = 2;
    opts.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL";
    Hamiltonian src_ham(opts);
    DenseHamiltonian dh(src_ham);
    const auto frm_dim = ci_utils::fermion_dim(src_ham.m_bd.m_nsite, src_ham.nelec());
    const auto bos_dim = ci_utils::boson_dim(src_ham.m_bd.m_nsite, src_ham.m_nboson_max, false);
    ASSERT_EQ(dh.m_ncol, frm_dim * bos_dim);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -0.38125085248276913);
}

TEST(DenseHamiltonian, BosonCouplingGeneralMaxOcc3) {
    fciqmc_config::Hamiltonian opts(nullptr);
    opts.m_fermion.m_fcidump.m_path = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    opts.m_ladder.m_ebdump.m_path = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL";
    opts.m_ladder.m_nboson_max = 3;
    opts.m_boson.m_bosdump.m_path = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL";
    Hamiltonian src_ham(opts);
    DenseHamiltonian dh(src_ham);
    const auto frm_dim = ci_utils::fermion_dim(src_ham.m_bd.m_nsite, src_ham.nelec());
    const auto bos_dim = ci_utils::boson_dim(src_ham.m_bd.m_nsite, src_ham.m_nboson_max, false);
    ASSERT_EQ(dh.m_ncol, frm_dim * bos_dim);
    auto solver = dh.diagonalize();
    ASSERT_FLOAT_EQ(solver.m_evals[0], -0.9998830020871416);
}
#endif