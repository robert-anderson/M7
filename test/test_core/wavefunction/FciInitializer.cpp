//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/FciInitializer.h>
#include <M7_lib/linalg/DenseHamiltonian.h>
#include <M7_lib/hamiltonian/frm/J1J2FrmHam.h>
#include "test_core/defs.h"

#ifdef ENABLE_FERMIONS
TEST(FciInitializer, N2) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    Hamiltonian ham(&frm_ham);
    FciInitOptions opt;
    opt.m_ritz_tol = 1e-7;
    DenseHamiltonian hmat(ham);
    v_t<ham_comp_t> dense_evals;
    dense::diag(hmat, dense_evals);
    ham_comp_t eval;
    auto results = FciInitializer::solve(ham, opt);
    results.get_eval(0, eval);
    ASSERT_NEAR_EQ(eval, dense_evals[0]);
}

TEST(FciInitializer, J1J2) {
    J1J2FrmHam frm_ham(0.25, lattice::make("ortho", {16}, {1}));
    Hamiltonian ham(&frm_ham);
    ham_comp_t eval;
    auto results = FciInitializer::solve(ham);
    results.get_eval(0, eval);
    ASSERT_NEAR_EQ(eval, -6.44708);
}
#endif

#ifdef ENABLE_BOSONS
TEST(FciInitializer, BosHub) {
    HubbardBosHam bos_ham(-0.1, lattice::make("ortho", {10}, {1}));
    Hamiltonian ham(&bos_ham);
    FciInitOptions opt;
    opt.m_nroot = 12;
    opt.m_ritz_tol = 1e-10;
    opt.m_diag_shift = -91.0;
    FciInitializer init(ham, opt);
}

TEST(FciInitializer, BosHubLoop) {
    for (ham_t u=-0.17; u<-0.14; u+=0.0005) {
        HubbardBosHam bos_ham(u, lattice::make("ortho", {10}, {1}));
        Hamiltonian ham(&bos_ham);
        FciInitOptions opt;
        opt.m_nroot = 10ul;
        opt.m_ritz_tol = 1e-10;
        opt.m_diag_shift = -91.0;
        FciInitializer init(ham, opt);
        auto results = init.get_results();
        results.bcast();
        hdf5::FileWriter fw(logging::format("bos_hub_u={:.4f}.h5", u));
        v_t<double> evals;
        results.get_evals(evals);
        fw.write_data("evals", evals);
        const ham_comp_t* evec = nullptr;
        results.get_evec(0ul, evec);
        const uintv_t shape = {results.nelement_evec(), results.nroot()};
        fw.write_data("evecs", evec, shape);
    }
}
#endif