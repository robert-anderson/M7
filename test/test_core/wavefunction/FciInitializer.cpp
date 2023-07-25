//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/CiInitializer.h>
#include <M7_lib/linalg/DenseHamiltonian.h>
#include <M7_lib/hamiltonian/frm/J1J2FrmHam.h>
#include "test_core/defs.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/mae/MaeTable.h"

#ifdef ENABLE_FERMIONS
#ifdef ENABLE_ARPACK
TEST(FciInitializer, N2Conns) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    Hamiltonian ham(&frm_ham);
    ci_init::Options opt;
    opt.m_ritz_tol = 1e-7;
    DenseHamiltonian hmat(ham);
    v_t<ham_comp_t> dense_evals;
    dense::diag(hmat, dense_evals);
    ham_comp_t eval;
    ci_init::FciSubspace subspace(&ham);
    auto results = ci_init::Initializer::solve(subspace, opt);
    results.get_eval(0, eval);
    ASSERT_NEAR_EQ(eval, dense_evals[0]);
}

TEST(FciInitializer, N2MbfPairs) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    Hamiltonian ham(&frm_ham);
    ci_init::Options opt;
    opt.m_ritz_tol = 1e-7;
    opt.m_loop_kind = ci_init::Options::MbfPairs;
    const ham_comp_t bench = -108.916561245698;
    ham_comp_t eval;
    ci_init::FciSubspace subspace(&ham);
    auto results = ci_init::Initializer::solve(subspace, opt);
    results.get_eval(0, eval);
    ASSERT_NEAR_EQ(eval, bench);
}

TEST(FciInitializer, N2Cisd) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    Hamiltonian ham(&frm_ham);
    ci_init::Options opt;
    opt.m_ritz_tol = 1e-7;
    opt.m_loop_kind = ci_init::Options::MbfPairs;
    const ham_comp_t bench = -108.89886691594;
    ham_comp_t eval;
    buffered::Mbf ref(ham.m_basis);
    mbf::set_aufbau_mbf(ref, ham.default_particles());
    ci_init::RefConnSubspace subspace(&ham, ref);
    auto results = ci_init::Initializer::solve(subspace, opt);
    results.get_eval(0, eval);
    ASSERT_NEAR_EQ(eval, bench);
}

TEST(FciInitializer, J1J2) {
    J1J2FrmHam frm_ham(0.25, lattice::make("ortho", {16}, {1}));
    Hamiltonian ham(&frm_ham);
    ham_comp_t eval;
    ci_init::FciSubspace subspace(&ham);
    auto results = ci_init::Initializer::solve(subspace);
    results.get_eval(0, eval);
    ASSERT_NEAR_EQ(eval, -6.44708);
}

#endif
#endif

#ifdef ENABLE_BOSONS
TEST(CiInitializer, BosHub) {
    HubbardBosHam bos_ham(-0.1, lattice::make("ortho", {10}, {1}));
    Hamiltonian ham(&bos_ham);
    CiInitOptions opt;
    opt.m_nroot = 12;
    opt.m_ritz_tol = 1e-10;
    opt.m_diag_shift = -91.0;
    CiInitializer init(ham, opt);
}

#if 0
TEST(CiInitializer, BosHubLoop) {
    const uint_t nsite = 9;
    const uint_t nbos = 9;
    {
        HubbardBosHam bos_ham(1.0, lattice::make("ortho", {nsite}, {1}));
        Hamiltonian ham(&bos_ham);
        const sys::Particles particles = {sys::frm::Electrons(0), sys::bos::Bosons(nbos, true)};
//        DenseHamiltonian dham(ham, particles);
//
//        hdf5::FileWriter fw1(logging::format("dense_{}site_{}bos.h5", nsite, nbos));
//        dham.save("data", fw1);

        auto iters = FciIters::make(ham, particles, false);
        buffered::BosOnv onv(ham.m_basis);
        v_t<buf_t> perms;
        perms.reserve(iters.m_single->m_niter*nsite);
        auto fn = [&onv, &perms]() {
            perms.insert(perms.end(), onv.cbegin(), onv.cend());
        };
        iters.m_single->loop(onv, fn);
        hdf5::FileWriter fw(logging::format("perms_{}site_{}bos.h5", nsite, nbos));
        fw.write_data("data", perms, {iters.m_single->m_niter, nsite});
    }

    for (ham_comp_t u=-0.5; u<0.2; u+=0.005) {
        HubbardBosHam bos_ham(u, lattice::make("ortho", {3, 3}, {1, 1}));
        Hamiltonian ham(&bos_ham);
        const sys::Particles particles = {sys::frm::Electrons(0), sys::bos::Bosons(nbos, true)};
        CiInitOptions opt;
        opt.m_nroot = nsite;
        opt.m_niter_max = 1000;
        opt.m_ritz_tol = 1e-8;
        opt.m_diag_shift = -20;
        auto results = CiInitializer::solve(ham, particles, opt);
        hdf5::FileWriter fw(logging::format("bos_hub_u={:.4f}.h5", u));
        v_t<double> evals;
        results.get_evals(evals);
        fw.write_data("evals", evals);
        uintv_t shape;
        shape.push_back(results.nelement_evec());//, results.nroot()};
        fw.write_data("evecs", results.get_evec(0ul), shape);
    }
}
#endif
#endif