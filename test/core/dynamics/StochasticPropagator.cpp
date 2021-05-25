//
// Created by RJA on 20/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/dynamics/Wavefunction.h"
#include "src/core/dynamics/StochasticPropagator.h"
#include "src/core/dynamics/Solver.h"


#ifndef ENABLE_BOSONS
TEST(StochasticPropagator, Test) {
    Options opts;
    opts.nwalker_initial = 100;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.05;
    opts.nwalker_target = 10000;
    opts.rdm_rank = 2;
    opts.replicate = true;
    opts.write_hdf5_fname = "rdm.h5";
    opts.ncycle_accumulate_mevs = 8000;
    opts.ncycle_mev_period = 10;
    opts.min_spawn_mag = 0.0;
    opts.min_death_mag = 0.2;
    opts.consolidate_spawns = false;
    opts.explicit_hf_conn_mevs = true;
    opts.output_mevs_periodically = true;
    opts.do_semistochastic = true;
    opts.ncycle_wait_mevs = 2200;
    opts.ncycle_wait_detsub = 2400;
    opts.init();

    FermionHamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv ref_onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        ref_onv.set({0, i});
        ref_onv.set({1, i});
    }

    Wavefunction wf(opts, ham.nsite());
    StochasticPropagator prop(ham, opts, wf.m_format);
    ASSERT_EQ(&wf.m_store.m_row.m_onv, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));

    auto ref_energy = ham.get_energy(ref_onv);
    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(opts.nwalker_initial);

    prop.m_shift.m_values = ref_energy+opts.shift_initial;
    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
}

TEST(StochasticPropagator, RdmTest) {
    Options opts;
    opts.nwalker_initial = 300.0;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.01;
    opts.load_balance_period = 5;
    opts.nload_balance_block_per_rank = 40;
    opts.nwalker_target = 200000;
    opts.shift_damp = 0.5;
    opts.shift_initial = 0.0;
    opts.ncycle_wait_mevs = 200;
    opts.ncycle_accumulate_mevs = 1000;
    opts.rdm_rank = 1;
    opts.replicate = true;
    opts.consolidate_spawns = true;
    opts.write_hdf5_fname = "test_rdm_save.h5";
    opts.init();

    //const auto benchmark = -99.9421389039331
    FermionHamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv ref_onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        ref_onv.set({0, i});
        ref_onv.set({1, i});
    }

    Wavefunction wf(opts, ham.nsite());
    ASSERT_EQ(wf.npart(), 2);
    StochasticPropagator prop(ham, opts, wf.m_format);
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ASSERT_EQ(&wf.m_store.m_row.m_onv, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));

    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, opts.nwalker_initial);
    wf.set_weight(1, opts.nwalker_initial);

    prop.m_shift.m_values = ref_energy+opts.shift_initial;
    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
}

TEST(StochasticPropagator, Hdf5) {
    Options opts;
    opts.nwalker_initial = 10.0;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.01;
    opts.load_balance_period = 5;
    opts.nload_balance_block_per_rank = 40;
    opts.nwalker_target = 50000;
    opts.shift_damp = 0.5;
    opts.shift_initial = 0.0;
    opts.ncycle = 50000;
    //opts.write_hdf5_fname = "test_wf_save.h5";
    opts.read_hdf5_fname = "test_wf_save.h5";
    opts.init();
    //const auto benchmark = -108.916561245585;
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv ref_onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        ref_onv.set({0, i});
        ref_onv.set({1, i});
    }

    Wavefunction wf(opts, ham.nsite());
    StochasticPropagator prop(ham, opts, wf.npart());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ASSERT_EQ(&wf.m_store.m_row.m_onv, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));

    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy+opts.shift_initial;
    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
}

TEST(StochasticPropagator, Hubbard) {
    Options opts;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.001;
    opts.nwalker_target = 20;
    opts.nwalker_initial = 1;
    opts.min_spawn_mag = 0.0;
    opts.min_death_mag = 0.0;
    opts.shift_damp = 1;
    opts.shift_update_period = 10;
    opts.ncycle = 2000000;
    opts.spf_uniform_twf = 1;
    opts.rdm_rank = 0;
    opts.init();

    // -4.22963352
    // -1.953145309
    Hamiltonian<> ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);

    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> ref_onv(ham.nsite());
    bool spin = false;
    //for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
    for (size_t i = 0ul; i < ham.nelec(); ++i) {
//        ref_onv.set({0, i});
//        ref_onv.set({1, i});
        ref_onv.set({spin, i});
        spin=!spin;
    }

    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    StochasticPropagator prop(ham, opts, wf.npart());
    auto ref_energy = ham.get_energy(ref_onv);
    prop.m_shift.m_values = ref_energy;

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, opts.nwalker_initial);

    Solver solver(prop, wf, ref_loc);
    std::cout << "Reference Energy: " << ref_energy << std::endl;
    solver.execute(opts.ncycle);
}

#else

TEST(StochasticPropagator, BosonTest) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 10000;
    opts.shift_damp = 0.4;
    opts.ncycle = 5000;
    opts.rdm_rank = 0;
    opts.replicate = false;
    //opts.spf_uniform_twf = true;
    //opts.write_hdf5_fname = "test_wf_save.h5";
    //opts.read_hdf5_fname = "test_wf_save.h5";
    opts.init();

    // -10.328242246088791
    Hamiltonian<> ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 2, 1.4, 0.3);

    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> ref_onv(ham.nsite());
    ham.set_hf_onv(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    StochasticPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);
    prop.m_shift.m_values = ref_energy;

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);
    Solver solver(prop, wf, ref_loc);

    solver.execute(opts.ncycle);
}
#endif
