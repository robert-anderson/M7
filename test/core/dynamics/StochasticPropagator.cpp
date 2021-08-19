//
// Created by RJA on 20/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/wavefunction/Wavefunction.h"
#include "src/core/dynamics/StochasticPropagator.h"
#include "src/core/dynamics/Solver.h"


#ifndef ENABLE_BOSONS
TEST(StochasticPropagator, Test) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 100;
    opts.m_propagator.m_nadd = 3.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 100000;
    opts.m_av_ests.m_rdm.m_ranks = {"2"};
    opts.m_av_ests.m_ncycle = 2000;
    opts.m_av_ests.m_stats_period = 10;
    opts.m_propagator.m_min_spawn_mag = 0.2;
    opts.m_propagator.m_min_death_mag = 0.2;
    opts.m_av_ests.m_delay = 1000;
    //opts.ncycle_wait_detsub = 10;
    opts.verify();

    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::FrmOnv ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());
    StochasticPropagator prop(ham, opts, wf.m_format);
    ASSERT_EQ(&wf.m_store.m_row.m_mbf, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));

    auto ref_energy = ham.get_energy(ref_onv);
    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(defs::wf_t(opts.m_wavefunction.m_nw_init));

    prop.m_shift.m_values = ref_energy+opts.m_shift.m_init;
    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}

TEST(StochasticPropagator, RdmTest) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 300.0;
    opts.m_propagator.m_nadd = 3.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_wavefunction.m_load_balancing.m_period = 5;
    opts.m_wavefunction.m_load_balancing.m_nblock_per_rank = 40;
    opts.m_propagator.m_nw_target = 200000;
    opts.m_shift.m_damp = 0.5;
    opts.m_av_ests.m_delay = 200;
    opts.m_av_ests.m_ncycle = 1000;
    opts.m_av_ests.m_rdm.m_ranks = {"1"};
    opts.verify();

    //const auto benchmark = -99.9421389039331
    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::FrmOnv ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());
    ASSERT_EQ(wf.npart(), 2);
    StochasticPropagator prop(ham, opts, wf.m_format);
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ASSERT_EQ(&wf.m_store.m_row.m_mbf, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));

    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, defs::wf_t(opts.m_wavefunction.m_nw_init));
    wf.set_weight(1, defs::wf_t(opts.m_wavefunction.m_nw_init));

    prop.m_shift.m_values = ref_energy+opts.m_shift.m_init;
    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}

TEST(StochasticPropagator, Hdf5) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10.0;
    opts.m_propagator.m_nadd = 3.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_wavefunction.m_load_balancing.m_period = 5;
    opts.m_wavefunction.m_load_balancing.m_nblock_per_rank = 40;
    opts.m_propagator.m_nw_target = 50000;
    opts.m_shift.m_damp = 0.5;
    opts.m_shift.m_init = 0.0;
    opts.m_propagator.m_ncycle = 50000;
    //opts.write_hdf5_fname = "test_wf_save.h5";
    opts.verify();
    //const auto benchmark = -108.916561245585;
    Hamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::FrmOnv ref_onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        ref_onv.set({0, i});
        ref_onv.set({1, i});
    }

    Wavefunction wf(opts, ham.nsite());
    StochasticPropagator prop(ham, opts, wf.npart());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ASSERT_EQ(&wf.m_store.m_row.m_mbf, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));

    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, defs::wf_t(ref_energy));

    prop.m_shift.m_values = ref_energy+opts.m_shift.m_init;
    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}

TEST(StochasticPropagator, Hubbard) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 1;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 100;
    opts.m_propagator.m_min_spawn_mag = 0.0;
    opts.m_propagator.m_min_death_mag = 0.0;
    opts.m_shift.m_damp = 0.1;
    opts.m_shift.m_period = 1;
    opts.m_shift.m_ncycle_av = 5000;
    opts.m_propagator.m_ncycle = 2000000;
    opts.m_inst_ests.m_spf_uniform_twf = false;
    opts.m_shift.m_reweight.m_ncycle = 200;
    opts.m_shift.m_reweight.m_delay = 20000;
    opts.m_av_ests.m_rdm.m_ranks = {};
    opts.verify();

    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);

    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    StochasticPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);
    prop.m_shift.m_values = ref_energy;

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, defs::wf_t(opts.m_wavefunction.m_nw_init));

    Solver solver(opts, prop, wf, ref_loc);
    std::cout << "Reference Energy: " << ref_energy << std::endl;
    solver.execute(opts.m_propagator.m_ncycle);
}

TEST(StochasticPropagator, ExcitedStates) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nroot = 3;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 3.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 10000;
    opts.verify();
    // -99.9421389039332
    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());

    if (wf.nreplica() && opts.m_wavefunction.m_nroot==1){ASSERT_EQ(wf.npart(), 1);}
    StochasticPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    for (size_t ipart=0ul; ipart<wf.npart(); ++ipart)
        wf.set_weight(ipart, defs::wf_t(opts.m_wavefunction.m_nw_init));


    buffered::Mbf excit_onv(ham.nsite());
    excit_onv = ref_onv;
    excit_onv.excite(ham.nelec()/2-1, ham.nelec()/2);
    wf.create_row(0, excit_onv, ham.get_energy(excit_onv), false);
    ASSERT(wf.m_store.m_row.index()==1ul);
    wf.set_weight(1, defs::wf_t(opts.m_wavefunction.m_nw_init));


    excit_onv = ref_onv;
    excit_onv.excite(ham.nsite()+ham.nelec()/2-1, ham.nsite()+ham.nelec()/2);
    wf.create_row(0, excit_onv, ham.get_energy(excit_onv), false);
    ASSERT(wf.m_store.m_row.index()==2ul);
    wf.set_weight(2, defs::wf_t(opts.m_wavefunction.m_nw_init));


    prop.m_shift.m_values = ref_energy;

    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
    //std::cout << solver.mevs().m_fermion_rdm->get_energy(ham.m_frm)-prop.m_shift.m_values[0]<< std::endl;
}


#else

TEST(StochasticPropagator, BosonTest) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 100;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.001;
    opts.m_propagator.m_nw_target = 10000;
    opts.m_shift.m_period = 1;
    opts.m_wavefunction.m_replicate = false;
    opts.m_propagator.m_min_spawn_mag = 0.2;
    opts.m_propagator.m_min_death_mag = 0.2;
    opts.m_propagator.m_consolidate_spawns = false;
    opts.verify();

    // -10.328242246088791
    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 0, 1.4, 0.3);

    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::Mbf ref(ham.nsite());
    ham.set_hf_mbf(ref, 0);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    StochasticPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref);
    prop.m_shift.m_values = ref_energy;

    auto ref_loc = wf.create_row(0, ref, ref_energy, 1);
    wf.set_weight(0, ref_energy);
    Solver solver(opts, prop, wf, ref_loc);

    solver.execute(opts.m_propagator.m_ncycle);
}
#endif
