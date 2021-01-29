//
// Created by rja on 10/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/dynamics/Wavefunction.h"
#include "src/core/dynamics/ExactPropagator.h"
#include "src/core/dynamics/Solver.h"


void print(const Wavefunction& wf){
    for (size_t i=0ul; i<wf.m_store.m_hwm; ++i){
//        std::cout << wf.m_walkers.m_onv(i).to_string() << " ";
//        std::cout << wf.m_walkers.m_weight(i, 0, 0) << std::endl;
    }
}

#ifdef ENABLE_BOSONS
TEST(ExactPropagator, BosonTest) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 100000;
    opts.shift_damp = 0.4;
    // nboson_cutoff 1: -6.9875779675355165
    // nboson_cutoff 2: -10.328242246088791
    Hamiltonian<1> ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 1, 1.4, 0.3);


    ASSERT_TRUE(ham.spin_conserving());
    elements::Onv<1> onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        onv.m_fonv.set(0, i);
        onv.m_fonv.set(1, i);
    }
    Wavefunction wf(opts, ham.nsite());
    wf.expand(10, 800);
    ExactPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(onv);
    prop.m_shift = ref_energy;//benchmark;
    Solver solver(prop, wf, onv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < 20000; ++i) {
        solver.execute();
        std::cout << i << " " << wf.m_walkers.m_hwm << " " << std::sqrt(wf.square_norm()) << std::endl;
        for (size_t irow = 0ul; irow < wf.m_walkers.m_hwm; ++irow) {
//            std::cout
//                    << wf.m_walkers.m_onv(irow).to_string() << " "
//                    << wf.m_walkers.m_weight(irow, 0, 0) << " "
//                    << wf.m_walkers.m_flags.m_reference_connection(irow)
//                    << std::endl;
        }
//        std::cout << "\n";
    }
/*
    for (size_t i = 0ul; i < wf.m_walkers.m_hwm; ++i) {
        std::cout
                << wf.m_walkers.m_onv(i).to_string() << " "
                << wf.m_walkers.m_weight(i, 0, 0) << " "
                << wf.m_walkers.m_flags.m_reference_connection(i)
                << std::endl;
    }
    */
}
#endif

#ifndef ENABLE_BOSONS
TEST(ExactPropagator, Test) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.05;
    opts.nwalker_target = 10000;
    const size_t nsite = 6;
    //const auto benchmark = -108.81138657563143;
    Hamiltonian<> ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    elements::Onv<> ref_onv(ham.nsite());
    for (size_t i=0ul; i<ham.nelec()/2; ++i){ref_onv.set(0, i); ref_onv.set(1, i);}
    Wavefunction wf(opts, nsite);
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ExactPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(ref_onv);
    prop.m_shift = ref_energy;//benchmark;

    auto ref_loc = wf.create_walker(ref_onv, opts.nwalker_initial, ref_energy, 1);
    prop.m_shift = ref_energy;
    Solver solver(prop, wf, ref_loc);
    for (size_t i = 0ul; i < opts.ncycle; ++i) {
        solver.execute();
    }
}

TEST(ExactPropagator, Cr2Test) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 10000;
    //const auto benchmark = -108.81138657563143;
    FermionHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    elements::FermionOnv ref_onv(ham.nsite());
    for (size_t i=0ul; i<ham.nelec()/2; ++i){ref_onv.set(0, i); ref_onv.set(1, i);}
    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(ref_onv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    auto ref_loc = wf.create_walker(ref_onv, opts.nwalker_initial, ref_energy, 1);
    prop.m_shift = ref_energy;
    Solver solver(prop, wf, ref_loc);
    for (size_t i = 0ul; i < opts.ncycle; ++i) {
        solver.execute();
    }
}
#endif