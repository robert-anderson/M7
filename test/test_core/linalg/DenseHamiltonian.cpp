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
    GeneralFrmHam frm_ham({defs::assets_root + "/RHF_N2_6o6e/FCIDUMP"}, true);
    Hamiltonian ham(&frm_ham);
    DenseHamiltonian hmat(ham);
    std::vector<defs::ham_t> evals;
    dense::diag(hmat, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -108.916561245585, 1e-8));
}

TEST(DenseHamiltonian, HeisenbergFrmHam) {
    /*
     * https://doi.org/10.1016/0378-4363(78)90115-8
     */
    defs::inds nsites = {4, 6, 8, 10, 12};
    std::vector<defs::ham_comp_t> energies = {-2.0, -2.8027756375, -3.6510934085, -4.515446354, -5.387390917};
    for (size_t i=0ul; i<nsites.size(); ++i){
        auto nsite = nsites[i];
        auto energy = energies[i];
        HeisenbergFrmHam frm_ham(1.0, lattice::make({Lattice::Ortho, {nsite}, {1}}));
        Hamiltonian ham(&frm_ham);
        std::vector<double> evals;
        auto particles = ham.default_particles();
        DenseHamiltonian hmat(ham, particles);
        dense::diag(hmat, evals);
        ASSERT_TRUE(consts::nearly_equal(evals[0], energy, 1e-8));
    }
}

TEST(DenseHamiltonian, HF) {
    GeneralFrmHam frm_ham({defs::assets_root + "/HF_RDMs/FCIDUMP"}, true);
    Hamiltonian ham(&frm_ham);
    DenseHamiltonian hmat(ham);
    std::vector<double> evals;
    dense::diag(hmat, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -99.9421389039332, 1e-8));
}

TEST(DenseHamiltonian, PyscfX2cCheck) {
    GeneralFrmHam frm_ham({defs::assets_root + "/H2O_X2C/FCIDUMP"}, true);
    Hamiltonian ham(&frm_ham);
    DenseHamiltonian hmat(ham);
    std::vector<double> evals;
    dense::diag(hmat, evals);
    // compare the ground and first excited states to PySCF's values
    ASSERT_TRUE(consts::nearly_equal(evals[0], -76.08150945314577, 1e-10));
}

TEST(DenseHamiltonian, Hubbard3Site) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    Hamiltonian ham(&frm_ham);
    std::vector<double> evals;
    DenseHamiltonian hmat(ham, ham.default_particles(4));
    dense::diag(hmat, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], 2.0, 1e-10));
}

TEST(DenseHamiltonian, Hubbard4Site) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    Hamiltonian ham(&frm_ham);
    std::vector<double> evals;
    DenseHamiltonian hmat(ham);
    dense::diag(hmat, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -1.9531453086749293, 1e-10));
}

TEST(DenseHamiltonian, Hubbard6Site) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {6}, {0}}));
    Hamiltonian ham(&frm_ham);
    std::vector<double> evals;
    DenseHamiltonian hmat(ham);
    dense::diag(hmat, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -3.0925653194551845, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinNoCoupling) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 0.0, 0);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.0);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], 2.0, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinNoFrequencyOccCutoff2) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 2);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.0);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    const auto nmode = ham_src.m_basis.m_bos.m_nmode;
    const auto& frm_basis = ham_src.m_basis.m_frm;
    for (size_t n = 0ul; n < nmode; ++n) {
        for (size_t p = 0ul; p < frm_basis.m_nspinorb; ++p) {
            const auto psite = frm_basis.isite(p);
            for (size_t q = 0ul; q < frm_basis.m_nspinorb; ++q) {
                const auto qsite = frm_basis.isite(q);
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
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -7.699484522379835, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinNoFrequencyOccCutoff3) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 3);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.0);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -11.07271962268484, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinOccCutoff2) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 2);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.3);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -6.692966463435127, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinOccCutoff1) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 1);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.3);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -3.1699561178752873, 1e-10));
}

TEST(DenseHamiltonian, HubbardHolsteinOccCutoff3) {
    HubbardFrmHam frm_ham(4.0, lattice::make({Lattice::Ortho, {3}, {0}}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 3);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.3);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -9.423844225360671, 1e-10));
}

TEST(DenseHamiltonian, BosonCouplingGeneralOccCutoff1) {
    conf::Hamiltonian opts(nullptr);
    GeneralFrmHam frm_ham({defs::assets_root + "/Hubbard_U4_3site/FCIDUMP"}, true);
    GeneralLadderHam frmbos_ham({defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 1);
    GeneralBosHam bos_ham({defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL"}, frmbos_ham.m_basis.m_bos.m_occ_cutoff);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    auto particles = ham_src.default_particles(4);
    DenseHamiltonian ham(ham_src, particles);
    std::vector<double> evals;
    dense::diag(ham, evals);
    std::cout << evals[0] << std::endl;
    ASSERT_TRUE(consts::nearly_equal(evals[0], 0.5090148148366922, 1e-10));
}

TEST(DenseHamiltonian, BosonCouplingGeneralOccCutoff2) {
    conf::Hamiltonian opts(nullptr);
    GeneralFrmHam frm_ham({defs::assets_root + "/Hubbard_U4_3site/FCIDUMP"}, true);
    GeneralLadderHam frmbos_ham({defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 2);
    GeneralBosHam bos_ham({defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL"}, frmbos_ham.m_basis.m_bos.m_occ_cutoff);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    auto particles = ham_src.default_particles(4);
    DenseHamiltonian ham(ham_src, particles);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -0.38125085248276913, 1e-10));
}

TEST(DenseHamiltonian, BosonCouplingGeneralOccCutoff3) {
    conf::Hamiltonian opts(nullptr);
    GeneralFrmHam frm_ham({defs::assets_root + "/Hubbard_U4_3site/FCIDUMP"}, true);
    GeneralLadderHam frmbos_ham({defs::assets_root + "/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 3);
    GeneralBosHam bos_ham({defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_GENERAL"}, frmbos_ham.m_basis.m_bos.m_occ_cutoff);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src);
    std::vector<double> evals;
    dense::diag(ham, evals);
    ASSERT_TRUE(consts::nearly_equal(evals[0], -0.9998830020871416, 1e-10));
}