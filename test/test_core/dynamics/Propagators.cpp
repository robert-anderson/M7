//
// Created by Robert J. Anderson on 29/06/2021.
//

#include <M7_lib/dynamics/Solver.h>
#include "gtest/gtest.h"
#include "M7_lib/dynamics/Propagators.h"


// tests of propagators should be done in the verification suite
#if 0

TEST(Propagators, BasicTest) {
    conf::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_nadd = 3.0;
    opts.m_propagator.m_tau_init = 0.01;
    opts.m_propagator.m_nw_target = 100000;
    opts.m_wavefunction.m_load_balancing.m_period = 1;
    opts.m_wavefunction.m_load_balancing.m_nblock_per_rank = 5;
    opts.verify();
    //const auto benchmark = -99.9421389039331
    //Hamiltonian ham(defs::assets_root + "/RHF_N2_CCPVDZ/FCIDUMP", false);
    Hamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    props::Stoch prop(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    if (ref_loc.is_mine())
        wf.set_weight(0, ref_energy);

    prop.m_shift.m_values = ref_energy;
    Solver solver(opts, prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}

TEST(Propagators, RefExcitTest) {
    conf::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_stochastic = false;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.05;
    opts.m_propagator.m_nw_target = 100;
    opts.m_av_ests.m_stats_period = 10;
    opts.m_av_ests.m_ref_excits.m_max_exlvl = 2;
    opts.m_av_ests.m_delay = 100;
    opts.m_av_ests.m_ncycle = 100;
    opts.verify();
    //const auto benchmark = -99.9421389039331
    Hamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::Mbf ref_onv(ham.nsite());
    ham.set_hf_mbf(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    auto prop = props::get(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, defs::wf_t(opts.m_wavefunction.m_nw_init));

    prop->m_shift.m_values = ref_energy;

    Solver solver(opts, *prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}
#endif