//
// Created by rja on 31/05/2021.
//

#include <src/core/io/Options.h>
#include <src/core/wavefunction/Wavefunction.h>
#include <src/core/dynamics/ExactPropagator.h>
#include <src/core/dynamics/Solver.h>
#include "gtest/gtest.h"


#ifndef ENABLE_BOSONS
/**
 *
 * @param opts
 *  program options
 * @return
 *  2RDM energy estimate
 */
defs::ham_comp_t fermion_rdm_energy_test(const fciqmc_config::Document& opts, bool explicit_hf_conns){
    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts, wf.m_format, explicit_hf_conns);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    for (size_t ipart = 0ul; ipart < wf.npart(); ++ipart)
        wf.set_weight(ipart, defs::wf_t(opts.m_wavefunction.m_nw_init));

    prop.m_shift.m_values = ref_energy+opts.m_shift.m_init.get();

    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
    return 0.0;//solver.mevs().m_fermion_rdm->get_energy(ham.m_frm);
}


void fermion_rdm_energy_opts(fciqmc_config::Document& opts){
    opts.m_archive.m_save_path = "M7.save";
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 50;
    opts.m_propagator.m_ncycle = 10000;
    opts.m_shift.m_init = 0.1;
    // enough cycles to solve for the ground state of this small system exactly
    opts.m_av_ests.m_delay = 100;
    opts.m_av_ests.m_ncycle = 400;
    opts.m_av_ests.m_stats_period = 10;
    opts.m_propagator.m_consolidate_spawns = false;
    opts.m_av_ests.m_fermion_rdm.m_ranks = {1, 2};
    opts.m_wavefunction.m_replicate = false;
}

TEST(FermionRdm, EnergyExactAverageEveryCycle) {
    fciqmc_config::Document opts;
    fermion_rdm_energy_opts(opts);
    // accumulate the average coefficient values every cycle
    opts.m_av_ests.m_stats_period = 1;
    opts.verify();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactAverageDivisibleNCycles) {
    fciqmc_config::Document opts;
    fermion_rdm_energy_opts(opts);
    // accumulate the average coefficient values every 10 cycles (100%10=0 so no need for finalization)
    opts.verify();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactAllLeftoverIters) {
    fciqmc_config::Document opts;
    fermion_rdm_energy_opts(opts);
    opts.m_av_ests.m_stats_period = 300;
    opts.verify();
    // 300>100 cycles so all contribs are averaged at the final loop over occupied ONVs
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactLeftoverIters) {
    fciqmc_config::Document opts;
    fermion_rdm_energy_opts(opts);
    opts.m_av_ests.m_stats_period = 13;
    /*
     * 100%13 = 9 cycles left over in averaging block when the ncycle_accumulate_mevs limit is reached, so a finalizing
     * loop is required for averaging correctness. this will be the default case for the remainder of the tests.
     */
    opts.verify();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactReplicated) {
    fciqmc_config::Document opts;
    fermion_rdm_energy_opts(opts);
    opts.m_av_ests.m_stats_period = 13;
    opts.m_wavefunction.m_replicate = true;
    opts.verify();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactExplicitHfConnections) {
    fciqmc_config::Document opts;
    fermion_rdm_energy_opts(opts);
    opts.m_av_ests.m_stats_period = 13;
    opts.m_wavefunction.m_replicate = true;
    opts.verify();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, true), -99.9421389039332);
}

#endif