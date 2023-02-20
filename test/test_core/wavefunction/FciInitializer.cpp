//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/FciInitializer.h>
#include <M7_lib/linalg/DenseHamiltonian.h>
#include <M7_lib/hamiltonian/frm/J1J2FrmHam.h>
#include "test_core/defs.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/mae/MaeTable.h"

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

TEST(FciInitializer, C4) {
    GeneralFrmHam frm_ham({"FCIDUMP"});
    Hamiltonian ham(&frm_ham);
    FciInitOptions opt;
    opt.m_ritz_tol = 1e-7;
    opt.m_niter_max = 1000;
    opt.m_nroot = 2;
    opt.m_diag_shift = -4;
    ham_comp_t eval;
    FciInitializer init(ham, opt);
    auto results = init.solve();


    buffered::MappedTable<RdmRow> c1("c1", {opsig::c_1100, 1ul});
    buffered::MappedTable<RdmRow> c2("c2", {opsig::c_2200, 1ul});
    buffered::MappedTable<RdmRow> c3("c3", {opsig::c_3300, 1ul});
    buffered::MappedTable<RdmRow> c4("c4", {opsig::c_4400, 1ul});

    wf_t nw_tot = 0.0;
    v_t<wf_t> nw_sectors(5, 0.0);

    if (mpi::i_am_root()) {

        results.get_eval(0, eval);
        logging::info("FCI energy: {}", eval);

        buffered::Mbf hf(ham.m_basis);
        const auto elecs = ham.default_particles().m_frm;
        mbf::set_aufbau_mbf(hf, elecs);
        conn::Mbf conn(hf);

        c1.resize(integer::combinatorial(
                2 * hf.m_basis.m_nsite - elecs, 1ul) * integer::combinatorial(uint_t(elecs), 1ul));
        c2.resize(integer::combinatorial(
                2 * hf.m_basis.m_nsite - elecs, 2ul) * integer::combinatorial(uint_t(elecs), 2ul));
        c3.resize(integer::combinatorial(
                2 * hf.m_basis.m_nsite - elecs, 3ul) * integer::combinatorial(uint_t(elecs), 3ul));
        c4.resize(integer::combinatorial(
                2 * hf.m_basis.m_nsite - elecs, 4ul) * integer::combinatorial(uint_t(elecs), 4ul));
        buffered::RdmInds c1_inds(opsig::c_1100);
        buffered::RdmInds c2_inds(opsig::c_2200);
        buffered::RdmInds c3_inds(opsig::c_3300);
        buffered::RdmInds c4_inds(opsig::c_4400);
        auto evec = results.get_evec(0);
        wf_t c0 = *evec;
        nw_sectors[0] = 1.0;
        const auto thresh = 1e-12;
        auto row = init.m_mbf_order_table.m_row;
        --evec;

        for (row.restart(); row; ++row) {
            ++evec;
            if (std::abs(c0) < std::abs(*evec)) {
                c0 = *evec;
                hf = row.m_field;
            }
        }
        logging::info("ref det: {}", hf.to_string());
        evec = results.get_evec(0);
        --evec;

        v_t<std::pair<RdmIndsField*, MappedTable<RdmRow>*>> pairs = {
                {nullptr,  nullptr},
                {&c1_inds, &c1},
                {&c2_inds, &c2},
                {&c3_inds, &c3},
                {&c4_inds, &c4}
        };

        for (row.restart(); row; ++row) {
            ++evec;
            if (std::abs(*evec) < thresh) continue;
            nw_tot += std::abs(*evec/c0);
            conn.connect(hf, row.m_field);
            const auto lvl = conn.exsig().nfrm_cre();
            if (!lvl) continue;
            if (lvl >= pairs.size()) continue;

            nw_sectors[lvl] += std::abs(*evec/c0);
            auto& inds = *pairs[lvl].first;
            auto& table = *pairs[lvl].second;
            inds = conn;
            table.insert(inds).m_values = *evec / c0;
        }
    }

    hdf5::FileWriter fw("ci.h5");
    c2.save(fw, "c2", mpi::i_am_root());
    c4.save(fw, "c4", mpi::i_am_root());
    hdf5::DatasetSaver::save_scalar(fw, "nw_tot", nw_tot);
    hdf5::DatasetSaver::save_vector(fw, "nw_sectors", nw_sectors);
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
            perms.insert(perms.end(), onv.begin(), onv.end());
        };
        iters.m_single->loop(onv, fn);
        hdf5::FileWriter fw(logging::format("perms_{}site_{}bos.h5", nsite, nbos));
        fw.write_data("data", perms, {iters.m_single->m_niter, nsite});
    }

    for (ham_comp_t u=-0.5; u<0.2; u+=0.005) {
        HubbardBosHam bos_ham(u, lattice::make("ortho", {3, 3}, {1, 1}));
        Hamiltonian ham(&bos_ham);
        const sys::Particles particles = {sys::frm::Electrons(0), sys::bos::Bosons(nbos, true)};
        FciInitOptions opt;
        opt.m_nroot = nsite;
        opt.m_niter_max = 1000;
        opt.m_ritz_tol = 1e-8;
        opt.m_diag_shift = -20;
        auto results = FciInitializer::solve(ham, particles, opt);
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