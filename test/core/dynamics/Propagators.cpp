//
// Created by rja on 29/06/2021.
//

#include <src/core/dynamics/Solver.h>
#include "gtest/gtest.h"
#include "src/core/dynamics/Propagators.h"

TEST(Propagators, AvCoeffTest) {
    fciqmc_config::Document opts;
    opts.m_wavefunction.m_nw_init = 10;
    opts.m_propagator.m_exact = true;
    opts.m_propagator.m_nadd = 0.0;
    opts.m_propagator.m_tau_init = 0.05;
    opts.m_propagator.m_nw_target = 100;
    opts.m_wavefunction.m_replicate = false;
    opts.m_observables.m_av_coeffs.m_max_exlvl = 2;
    opts.m_observables.m_delay = 100;
    opts.m_observables.m_ncycle = 100;
    opts.verify();
    //const auto benchmark = -99.9421389039331
    Hamiltonian<> ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> ref_onv(ham.nsite());
    ham.set_hf_onv(ref_onv, 0);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    auto prop = props::get(ham, opts, wf.m_format);
    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(0, opts.m_wavefunction.m_nw_init);

    prop->m_shift.m_values = ref_energy;

    Solver solver(opts, *prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}