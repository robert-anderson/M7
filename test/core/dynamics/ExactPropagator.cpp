//
// Created by rja on 10/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/dynamics/Wavefunction.h"
#include "src/core/dynamics/ExactPropagator.h"
#include "src/core/dynamics/Solver.h"


void print(const Wavefunction& wf){
    for (size_t i=0ul; i<wf.m_walkers.m_hwm; ++i){
        std::cout << wf.m_walkers.m_onv(i).to_string() << " ";
        std::cout << wf.m_walkers.m_weight(i, 0, 0) << std::endl;
    }
}

TEST(ExactPropagator, Test) {
    Options opts;
    opts.nwalker_initial = 100;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.05;
    fields::FermiBosOnv::params_t params{6, 6};
    const auto benchmark = -108.81138657563143;
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    elements::FermionOnv fonv(ham.nsite());
    for (size_t i=0ul; i<ham.nelec()/2; ++i){fonv.set(0, i); fonv.set(1, i);}
    Wavefunction wf(opts, params);
    wf.expand(1000, 8000);
    ExactPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(fonv);
    prop.m_shift = benchmark;
    Solver solver(prop, wf, fonv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    print(wf);

    for (size_t i=0ul; i<10000; ++i){
        solver.execute();
        std::cout << i << " " << wf.m_walkers.m_hwm << " " << std::sqrt(wf.square_norm()) << " " << benchmark-wf.energy(ham) << std::endl;
    }


}