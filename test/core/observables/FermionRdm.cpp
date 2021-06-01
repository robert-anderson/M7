//
// Created by rja on 31/05/2021.
//

#include <src/core/io/Options.h>
#include <src/core/dynamics/Wavefunction.h>
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
defs::ham_comp_t fermion_rdm_energy_test(const Options& opts){
    Hamiltonian<> ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    buffered::Onv<> ref_onv(ham.nsite());
    ham.set_hf_onv(ref_onv, 0);

    Wavefunction wf(opts, ham.nsite());
    ExactPropagator prop(ham, opts, wf.m_format, opts.explicit_hf_conn_mevs);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    for (size_t ipart = 0ul; ipart < wf.npart(); ++ipart) wf.set_weight(ipart, opts.nwalker_initial);

    prop.m_shift.m_values = ref_energy;

    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
    return solver.mevs().m_fermion_rdm->get_energy(ham);
}


void fermion_rdm_energy_opts(Options& opts){
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 0.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 50;
    // enough cycles to solve for the ground state of this small system exactly
    opts.ncycle_wait_mevs = 2000;
    opts.ncycle_accumulate_mevs = 100;
    opts.ncycle_mev_period = 10;
    opts.consolidate_spawns = false;
    opts.explicit_hf_conn_mevs = false;
    opts.rdm_rank = 2;
    opts.replicate = false;
}

TEST(FermionRdm, EnergyExactAverageEveryCycle) {
    Options opts;
    fermion_rdm_energy_opts(opts);
    // accumulate the average coefficient values every cycle
    opts.ncycle_mev_period = 1;
    opts.init();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactAverageDivisibleNCycles) {
    Options opts;
    fermion_rdm_energy_opts(opts);
    // accumulate the average coefficient values every 10 cycles (100%10=0 so no need for finalization)
    opts.init();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactAllLeftoverIters) {
    Options opts;
    fermion_rdm_energy_opts(opts);
    opts.ncycle_mev_period = 300;
    opts.init();
    // 300>100 cycles so all contribs are averaged at the final loop over occupied ONVs
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactLeftoverIters) {
    Options opts;
    fermion_rdm_energy_opts(opts);
    opts.ncycle_mev_period = 13;
    /*
     * 100%13 = 9 cycles left over in averaging block when the ncycle_accumulate_mevs limit is reached, so a finalizing
     * loop is required for averaging correctness. this will be the default case for the remainder of the tests.
     */
    opts.init();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactReplicated) {
    Options opts;
    fermion_rdm_energy_opts(opts);
    opts.ncycle_mev_period = 13;
    opts.replicate = true;
    opts.init();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactExplicitHfConnections) {
    Options opts;
    fermion_rdm_energy_opts(opts);
    opts.ncycle_mev_period = 13;
    opts.replicate = true;
    opts.explicit_hf_conn_mevs = true;
    opts.init();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts), -99.9421389039332);
}

#endif