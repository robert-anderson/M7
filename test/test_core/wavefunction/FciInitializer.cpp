//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/FciInitializer.h>
#include <M7_lib/linalg/DenseHamiltonian.h>
#include <M7_lib/hamiltonian/frm/J1J2FrmHam.h>
#include "test_core/defs.h"

#ifdef ENABLE_FERMIONS
TEST(FciInitializer, N2) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP", true});
    Hamiltonian ham(&frm_ham);
    ArnoldiOptions opt;
    opt.m_ritz_tol = 1e-7;
    FciInitializer init(ham, 0.0, opt);
    DenseHamiltonian hmat(ham);
    v_t<ham_t> dense_evals;
    dense::diag(hmat, dense_evals);
    ham_comp_t eval;
    auto results = init.get_results();
    results.bcast();
    results.get_eval(0, eval);
    ASSERT_NEARLY_EQ(eval, dense_evals[0]);
}

TEST(FciInitializer, J1J2) {
    J1J2FrmHam frm_ham(0.25, lattice::make("ortho", {16}, {1}));
    Hamiltonian ham(&frm_ham);
    FciInitializer init(ham, 0.0);
    ham_comp_t eval;
    auto results = init.get_results();
    results.bcast();
    results.get_eval(0, eval);
    ASSERT_NEARLY_EQ(eval, -6.44708);
}
#endif

#ifdef ENABLE_BOSONS
TEST(FciInitializer, BosHub) {
    HubbardBosHam bos_ham(-0.1, lattice::make("ortho", {10}, {1}));
    Hamiltonian ham(&bos_ham);
    ArnoldiOptions opt;
    opt.m_nroot = 12;
    opt.m_ritz_tol = 1e-10;
    FciInitializer init(ham, -91, opt);
}
#endif