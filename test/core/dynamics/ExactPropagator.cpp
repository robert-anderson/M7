//
// Created by rja on 10/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/dynamics/Wavefunction.h"
#include "src/core/dynamics/ExactPropagator.h"
#include "src/core/dynamics/Solver.h"



#if 0 
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
        onv.m_frm.set(0, i);
        onv.m_frm.set(1, i);
    }
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ExactPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(onv);


    auto ref_loc = wf.create_walker(onv, opts.nwalker_initial, ref_energy, 1);
    prop.m_values = ref_energy;
    Solver solver(prop, wf, ref_loc);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < 20000; ++i) {
        solver.execute();
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
    opts.rdm_rank = 1;
    opts.replicate = false;
    opts.write_hdf5_fname = "rdm.h5";
    opts.ncycle = 10000;
    //const auto benchmark = -108.916561245585
    Hamiltonian<> ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> ref_onv(ham.nsite());
    ham.set_hf_onv(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;

    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
}

TEST(ExactPropagator, RdmTest) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.05;
    opts.nwalker_target = 10000;
    //opts.rdm_rank = 2;
    opts.replicate = false;
    //const auto benchmark = -99.9421389039331
    Hamiltonian<> ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> ref_onv(ham.nsite());
    ham.set_hf_onv(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ExactPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;

    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
}



TEST(ExactPropagator, Hubbard) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 10000;
    opts.shift_damp = 0.4;
    opts.ncycle = 3000;
    opts.spf_uniform_twf = 0;
    opts.rdm_rank = 0;
    opts.init();

    Hamiltonian<> ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);

    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> ref_onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        ref_onv.set({0, i});
        ref_onv.set({1, i});
    }
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ExactPropagator prop(ham, opts, wf.npart());
    auto ref_energy = ham.get_energy(ref_onv);
    prop.m_shift.m_values = ref_energy;

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    Solver solver(prop, wf, ref_loc);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < 20000; ++i) {
        solver.execute();
    }
}

TEST(ExactPropagator, Cr2Test) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 10000;
    opts.ncycle = 50;
    //const auto benchmark = -108.916561245585;
    FermionHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv ref_onv(ham.nsite());
    for (size_t i=0ul; i<ham.nelec()/2; ++i){ref_onv.set(0, i); ref_onv.set({1, i});}
    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts, wf.npart());
    auto ref_energy = ham.get_energy(ref_onv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;
    Solver solver(prop, wf, ref_loc);

    std::cout <<
        wf.m_store.nbucket()
    << std::endl;

    for (size_t i = 0ul; i < opts.ncycle; ++i) {
        solver.execute();
        std::cout << wf.m_store.m_hwm << " " << wf.m_store.m_nrow << std::endl;
    }
}
#endif