//
// Created by Robert John Anderson on 2020-01-24.
//

#include "test_core/defs.h"
#include "M7_lib/linalg/DenseHamiltonian.h"

/*
 * exact diagonalization in the entire Hilbert space for integration testing of matrix elements for very small systems
 */

#ifdef ENABLE_COMPLEX
TEST(DenseHamiltonian, FciEnergyCheck4c) {
    DenseHamiltonian ham(Hamiltonian(PROJECT_ROOT"/assets/DHF_Be_STO-3G/FCIDUMP", false));
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(dtype::floats_nearly_equal(solver.m_evals[0], -14.40597603432, 1e-10));
    ASSERT_TRUE(dtype::floats_nearly_equal(solver.m_evals[1], -14.28883698406, 1e-10));
}
#endif
TEST(DenseHamiltonian, N2Rhf) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"}, true);
    Hamiltonian ham(&frm_ham);
    DenseHamiltonian hmat(ham);
    v_t<ham_t> evals;
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(evals[0], -108.916561245585);
}

TEST(DenseHamiltonian, HeisenbergFrmHam) {
    /*
     * https://doi.org/10.1016/0378-4363(78)90115-8
     */
    uintv_t nsites = {4, 6, 8, 10, 12};
    v_t<ham_comp_t> energies = {-2.0, -2.8027756375, -3.6510934085, -4.515446354, -5.387390917};
    for (uint_t i=0ul; i<nsites.size(); ++i){
        auto nsite = nsites[i];
        auto energy = energies[i];
        HeisenbergFrmHam frm_ham(1.0, lattice::make("ortho", {nsite}, {1}));
        Hamiltonian ham(&frm_ham);
        v_t<double> evals;
        auto particles = ham.default_particles();
        DenseHamiltonian hmat(ham, particles);
        dense::diag(hmat, evals);
        ASSERT_NEARLY_EQ(evals[0], energy);
    }
}

TEST(DenseHamiltonian, HF) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/HF_RDMs/FCIDUMP"}, true);
    Hamiltonian ham(&frm_ham);
    auto particles = ham.default_particles();
    ASSERT_EQ(uint_t(particles.m_frm), 6ul);
    DenseHamiltonian hmat(ham, particles);
    v_t<double> evals;
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(evals[0], -99.9421389039332);
}


TEST(DenseHamiltonian, N2Molcas) {
    /*
     * HF:      -108.9540866268
     * CASCI:   -109.02180323
     */
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5"}, true);
    Hamiltonian ham(&frm_ham);
    auto particles = ham.default_particles();
    ASSERT_EQ(uint_t(particles.m_frm), 6ul);
    DenseHamiltonian hmat(ham, particles);
    v_t<double> evals;
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(evals[0], -109.02180323);
}

TEST(DenseHamiltonian, PyscfX2cCheck) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/H2O_X2C/FCIDUMP"}, true);
    Hamiltonian ham(&frm_ham);
    DenseHamiltonian hmat(ham);
    v_t<double> evals;
    dense::diag(hmat, evals);
    // compare the ground and first excited states to PySCF's values
    ASSERT_NEARLY_EQ(evals[0], -76.08150945314577);
}

TEST(DenseHamiltonian, Hubbard3Site) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    Hamiltonian ham(&frm_ham);
    v_t<double> evals;
    DenseHamiltonian hmat(ham, ham.default_particles(4));
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(evals[0], 2.0);
}

TEST(DenseHamiltonian, Hubbard4Site) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {4}, {0}));
    Hamiltonian ham(&frm_ham);
    v_t<double> evals;
    DenseHamiltonian hmat(ham);
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(evals[0], -1.9531453086749293);
}

TEST(DenseHamiltonian, Hubbard6Site) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {6}, {0}));
    Hamiltonian ham(&frm_ham);
    v_t<double> evals;
    DenseHamiltonian hmat(ham);
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(evals[0], -3.0925653194551845);
}

TEST(DenseHamiltonian, HubbardHolsteinNoCoupling) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    ASSERT_EQ(frm_ham.m_basis.m_nsite, 3ul);
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 0.0, 0);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.0);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    ASSERT_EQ(ham_src.m_basis.m_frm.m_nsite, 3ul);
    auto particles = ham_src.default_particles(4);
    DenseHamiltonian ham(ham_src, particles);
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], 2.0);
}

TEST(DenseHamiltonian, HubbardHolsteinNoFrequencyOccCutoff2) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 2);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.0);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    const auto nmode = ham_src.m_basis.m_bos.m_nmode;
    const auto& frm_basis = ham_src.m_basis.m_frm;
    for (uint_t n = 0ul; n < nmode; ++n) {
        for (uint_t p = 0ul; p < frm_basis.m_nspinorb; ++p) {
            const auto psite = frm_basis.isite(p);
            for (uint_t q = 0ul; q < frm_basis.m_nspinorb; ++q) {
                const auto qsite = frm_basis.isite(q);
                if (n == psite && psite == qsite) {
                    ASSERT_NEARLY_EQ(arith::real(ham_src.m_frmbos.get_coeff_1101(n, p, q)), 1.4);
                    ASSERT_NEARLY_EQ(arith::real(ham_src.m_frmbos.get_coeff_1110(n, p, q)), 1.4);
                } else {
                    ASSERT_NEARLY_EQ(arith::real(ham_src.m_frmbos.get_coeff_1101(n, p, q)), 0.0);
                    ASSERT_NEARLY_EQ(arith::real(ham_src.m_frmbos.get_coeff_1110(n, p, q)), 0.0);
                }
            }
        }
    }
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -7.699484522379835);
}

TEST(DenseHamiltonian, HubbardHolsteinNoFrequencyOccCutoff3) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 3);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.0);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -11.07271962268484);
}

TEST(DenseHamiltonian, HubbardHolsteinOccCutoff2) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 2);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.3);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -6.692966463435127);
}

TEST(DenseHamiltonian, HubbardHolsteinOccCutoff1) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 1);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.3);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -3.1699561178752873);
}

TEST(DenseHamiltonian, HubbardHolsteinOccCutoff3) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {3}, {0}));
    HolsteinLadderHam frmbos_ham(frm_ham.m_basis, 1.4, 3);
    NumOpBosHam bos_ham(frmbos_ham.m_basis.m_bos, 0.3);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src, ham_src.default_particles(4));
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -9.423844225360671);
}

TEST(DenseHamiltonian, GeneralFrmBosOccCutoff1) {
    conf::Hamiltonian opts(nullptr);
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/FCIDUMP"}, true);
    GeneralLadderHam frmbos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 1);
    GeneralBosHam bos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/BOSDUMP_GENERAL"}, frmbos_ham.m_basis.m_bos.m_occ_cutoff);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    auto particles = ham_src.default_particles(4);
    DenseHamiltonian ham(ham_src, particles);
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], 0.5090148148366922);
}

TEST(DenseHamiltonian, GeneralFrmBosOccCutoff2) {
    conf::Hamiltonian opts(nullptr);
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/FCIDUMP"}, true);
    GeneralLadderHam frmbos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 2);
    GeneralBosHam bos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/BOSDUMP_GENERAL"}, frmbos_ham.m_basis.m_bos.m_occ_cutoff);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    auto particles = ham_src.default_particles(4);
    DenseHamiltonian ham(ham_src, particles);
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -0.38125085248276913);
}

TEST(DenseHamiltonian, GeneralFrmBosOccCutoff3) {
    conf::Hamiltonian opts(nullptr);
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/FCIDUMP"}, true);
    GeneralLadderHam frmbos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 3);
    GeneralBosHam bos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/BOSDUMP_GENERAL"}, frmbos_ham.m_basis.m_bos.m_occ_cutoff);
    Hamiltonian ham_src(&frm_ham, &frmbos_ham, &bos_ham);
    DenseHamiltonian ham(ham_src);
    v_t<double> evals;
    dense::diag(ham, evals);
    ASSERT_NEARLY_EQ(evals[0], -0.9998830020871416);
}