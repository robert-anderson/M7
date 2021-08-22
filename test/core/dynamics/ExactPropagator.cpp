//
// Created by rja on 10/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/wavefunction/Wavefunction.h"
#include "src/core/dynamics/ExactPropagator.h"
#include "src/core/dynamics/Solver.h"


#ifdef ENABLE_BOSONS
TEST(ExactPropagator, BosonTest) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_wavefunction.m_nw_init = 100;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 10000;
    opts.m_shift.m_period = 1;
    opts.m_wavefunction.m_replicate = false;
    opts.m_propagator.m_min_spawn_mag = 0.2;
    opts.m_propagator.m_min_death_mag = 0.2;
    opts.m_propagator.m_consolidate_spawns = false;
    opts.verify();
    // nboson_cutoff 1: -6.9875779675355165
    // nboson_cutoff 2: -10.328242246088791
    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 1, 1.4, 0.3);
    ASSERT_TRUE(ham.m_frm.spin_conserving());
    buffered::Mbf ref(ham.nsite());
    ham.set_hf_mbf(ref, 0);

    Wavefunction wf(opts, ham.nsite());
    if (!opts.m_wavefunction.m_replicate && opts.m_wavefunction.m_nroot == 1) { ASSERT_EQ(wf.npart(), 1); }
    bool explicit_hf_conn_mevs = false; //TODO
    ExactPropagator prop(ham, opts, wf.m_format, explicit_hf_conn_mevs);
    auto ref_energy = ham.get_energy(ref);

    auto ref_loc = wf.create_row(0, ref, ref_energy, 1);
    for (size_t ipart = 0ul; ipart < wf.npart(); ++ipart) wf.set_weight(ipart, opts.m_wavefunction.m_nw_init);

    prop.m_shift.m_values = ref_energy;

    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
    //std::cout << solver.mevs().m_fermion_rdm->get_energy(ham.m_frm) - prop.m_shift.m_values[0] << std::endl;
}

#else

TEST(ExactPropagator, DeterministicSubspace) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.05;
    opts.m_propagator.m_nw_target = 1000;
    opts.m_wavefunction.m_load_balancing.m_period = 1;
    opts.m_wavefunction.m_load_balancing.m_nblock_per_rank = 5;
    opts.verify();
    //const auto benchmark = -99.9421389039331
    //Hamiltonian<> ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    Hamiltonian ham(defs::assets_root + "/RHF_N2_CCPVDZ/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    if (ref_loc.is_mine())
        wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;
    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}


TEST(ExactPropagator, Test) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nroot = 3;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 50;
    opts.m_av_ests.m_delay = 4000;
    opts.m_av_ests.m_ncycle = 200;
    opts.m_av_ests.m_stats_period = 13;
    opts.m_av_ests.m_rdm.m_ranks = {};
    opts.verify();
    // -99.9421389039332
    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());

    if (wf.nreplica()==1 && opts.m_wavefunction.m_nroot==1){ASSERT_EQ(wf.npart(), 1);}
    ExactPropagator prop(ham, opts, wf.m_format, true);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    for (size_t ipart=0ul; ipart<wf.npart(); ++ipart) wf.set_weight(ipart, defs::wf_t(opts.m_wavefunction.m_nw_init));


    buffered::Mbf excit_onv(ham.nsite());
    excit_onv = ref_onv;
    excit_onv.excite(ham.nelec()/2-1, ham.nelec()/2);
    wf.create_row(0, excit_onv, ham.get_energy(excit_onv), false);
    ASSERT(wf.m_store.m_row.index()==1);
    wf.set_weight(1, defs::wf_t(opts.m_wavefunction.m_nw_init));


    excit_onv = ref_onv;
    excit_onv.excite(ham.nsite()+ham.nelec()/2-1, ham.nsite()+ham.nelec()/2);
    wf.create_row(0, excit_onv, ham.get_energy(excit_onv), false);
    ASSERT(wf.m_store.m_row.index()==2ul);
    wf.set_weight(2, defs::wf_t(opts.m_wavefunction.m_nw_init));

    std::cout <<
              wf.m_store.to_string()
    << std::endl;
    prop.m_shift.m_values = ref_energy;

    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
    //std::cout << solver.mevs().m_fermion_rdm->get_energy(ham.m_frm)-prop.m_shift.m_values[0]<< std::endl;
}

TEST(ExactPropagator, RdmTest) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.05;
    opts.m_propagator.m_nw_target = 1000;
    opts.verify();
    //const auto benchmark = -99.9421389039331
    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ExactPropagator prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;

    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}



TEST(ExactPropagator, Hubbard) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 10000;
    opts.m_shift.m_damp = 0.4;
    opts.m_propagator.m_ncycle = 3000;
    opts.m_inst_ests.m_spf_uniform_twf = 0;
    opts.m_av_ests.m_rdm.m_ranks = {};
    opts.verify();

    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);

    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::Mbf ref_onv(ham.nsite());
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

    Solver solver(opts, prop, wf, ref_loc);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < 20000; ++i) {
        solver.execute(opts.m_propagator.m_ncycle);
    }
}

TEST(ExactPropagator, Cr2Test) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_reference.m_redef_thresh = 2.0;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 10000;
    opts.m_propagator.m_ncycle = 50;
    opts.verify();
    //const auto benchmark = -108.916561245585;
    Hamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::FrmOnv ref_onv(ham.nsite());
    for (size_t i=0ul; i<ham.nelec()/2; ++i){ref_onv.set(0, i); ref_onv.set({1, i});}
    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts, wf.npart());
    auto ref_energy = ham.get_energy(ref_onv);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;
    Solver solver(opts, prop, wf, ref_loc);

    std::cout <<
        wf.m_store.nbucket()
    << std::endl;

    solver.execute(opts.m_propagator.m_ncycle);
}
#endif