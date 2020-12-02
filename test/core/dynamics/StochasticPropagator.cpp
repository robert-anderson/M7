//
// Created by RJA on 20/11/2020.
//


#include "gtest/gtest.h"
#include "src/core/dynamics/Wavefunction.h"
#include "src/core/dynamics/StochasticPropagator.h"
#include "src/core/dynamics/Solver.h"

/*
TEST(StochasticPropagator, Test) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.05;
    opts.nwalker_target = 100000;
//const auto benchmark = -108.81138657563143;
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    elements::FermionOnv fonv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        fonv.set(0, i);
        fonv.set(1, i);
    }
    Wavefunction wf(opts, ham.nsite());
    wf.expand(10, 800);
    StochasticPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(fonv);
    prop.m_shift = ref_energy;//benchmark;
    Solver solver(prop, wf, fonv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < 10000; ++i) {
        solver.execute();
        std::cout << i << " " << wf.m_walkers.m_hwm << " " << std::sqrt(wf.square_norm()) << std::endl;
    }
}
 */

TEST(StochasticPropagator, BosonTest) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 10000;
    opts.shift_damp = 0.4;
    opts.ncycle = 3000;


    // -10.328242246088791

    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 2, 1.4, 0.3);

    ASSERT_TRUE(ham.spin_conserving());
    elements::Onv onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        onv.m_fonv.set(0, i);
        onv.m_fonv.set(1, i);
    }
    Wavefunction wf(opts, ham.nsite());
    wf.expand(10, 800);
    StochasticPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(onv);
    prop.m_shift = ref_energy;//benchmark;
    Solver solver(prop, wf, onv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < opts.ncycle; ++i) {
        solver.execute();
        std::cout << i << " " << wf.m_walkers.m_hwm << " " << std::sqrt(wf.square_norm()) << std::endl;
    }

    for (size_t i = 0ul; i < wf.m_walkers.m_hwm; ++i) {
        std::cout
        << wf.m_walkers.m_onv(i).to_string() << " "
        << wf.m_walkers.m_weight(i, 0, 0)
        << std::endl;
    }
}