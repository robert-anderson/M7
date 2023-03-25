//
// Created by Robert J. Anderson on 31/05/2021.
//

#include <M7_lib/io/Options.h>
#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/propagator/ExactLinear.h>
#include <M7_lib/dynamics/Solver.h>
#include <M7_lib/field/Mbf.h>
#include "gtest/gtest.h"


#ifndef ENABLE_BOSONS
/**
 *
 * @param opts
 *  program options
 * @return
 *  2RDM energy estimate
 */
ham_comp_t fermion_rdm_energy_test(const conf::Document& opts, bool explicit_hf_conns){
    Hamiltonian ham({opts.m_hamiltonian, opts.m_basis, opts.m_particles});
    auto particles = ham.default_particles();
    buffered::Mbf ref_onv(ham.m_basis);
    mbf::set_aufbau_mbf(ref_onv, particles);

    wf::Vectors wf(opts, ham);
    ExactLinear prop(ham, opts, wf, explicit_hf_conns);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, 1);
    if (ref_loc.is_mine()) {
        auto walker = wf.m_store.m_row;
        walker.jump(ref_loc.m_irec);
        for (uint_t ipart = 0ul; ipart < wf.npart(); ++ipart)
            wf.set_weight(walker, ipart, wf_t(opts.m_wavefunction.m_nw_init));
    }

    prop.m_shift.m_values = ref_energy+opts.m_shift.m_init;

    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
    return solver.m_maes.m_rdms.get_energy(ham);
}


void fermion_rdm_energy_opts(conf::Document& opts){
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = PROJECT_ROOT"/assets/HF_RDMs/FCIDUMP";
    opts.m_av_ests.m_rdm.m_save.m_path = "M7.rdm.save.h5";
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
    opts.m_av_ests.m_rdm.m_ranks = {"1", "2"};
}

TEST(FermionRdm, EnergyExactAverageEveryCycle) {
    conf::Document opts;
    fermion_rdm_energy_opts(opts);
    // accumulate the average coefficient values every cycle
    opts.m_av_ests.m_stats_period = 1;
    opts.validate();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactAverageDivisibleNCycles) {
    conf::Document opts;
    fermion_rdm_energy_opts(opts);
    // accumulate the average coefficient values every 10 cycles (100%10=0 so no need for finalization)
    opts.validate();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactAllLeftoverIters) {
    conf::Document opts;
    fermion_rdm_energy_opts(opts);
    opts.m_av_ests.m_stats_period = 300;
    opts.validate();
    // 300>100 cycles so all contribs are averaged at the final loop over occupied ONVs
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

TEST(FermionRdm, EnergyExactLeftoverIters) {
    conf::Document opts;
    fermion_rdm_energy_opts(opts);
    opts.m_av_ests.m_stats_period = 13;
    /*
     * 100%13 = 9 cycles left over in averaging block when the ncycle_accumulate_mevs limit is reached, so a finalizing
     * loop is required for averaging correctness. this will be the default case for the remainder of the tests.
     */
    opts.validate();
    ASSERT_FLOAT_EQ(fermion_rdm_energy_test(opts, false), -99.9421389039332);
}

#endif