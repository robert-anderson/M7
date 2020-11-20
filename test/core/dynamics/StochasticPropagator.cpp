//
// Created by RJA on 20/11/2020.
//


#include "gtest/gtest.h"
#include "src/core/dynamics/Wavefunction.h"
#include "src/core/dynamics/StochasticPropagator.h"
#include "src/core/dynamics/Solver.h"

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

    for (size_t i = 0ul; i < 1000; ++i) {
        solver.execute();
        std::cout << i << " " << wf.m_walkers.m_hwm << " " << std::sqrt(wf.square_norm()) << std::endl;
//std::cout << solver.reference().proj_energy()-benchmark << std::endl;
    }
}