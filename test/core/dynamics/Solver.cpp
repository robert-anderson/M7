//
// Created by rja on 14/07/2021.
//

#include <src/core/dynamics/Propagators.h>
#include "gtest/gtest.h"
#include "src/core/dynamics/Solver.h"

TEST(Solver, RecvSort) {
    fciqmc_config::Document opts;
    opts.m_propagator.m_nw_target = 1000;
    opts.m_propagator.m_consolidate_spawns = true;
    opts.m_wavefunction.m_nw_init = 10;
    opts.verify();
    const size_t nsite = 4;
    Wavefunction wf(opts, nsite);

    /*
     * set up some arbitrary system of determinants
     */
    buffered::FrmOnv onv0(nsite);
    onv0 = {0, 5, 6, 7};
    buffered::FrmOnv onv1(nsite);
    onv1 = {1, 5, 6, 7};
    buffered::FrmOnv onv2(nsite);
    onv2 = {2, 5, 6, 7};
    buffered::FrmOnv onv3(nsite);
    onv3 = {3, 5, 6, 7};

    auto& recv = wf.m_comm.recv();
    auto& row = recv.m_row;
    row.push_back_jump();
    row.m_dst_onv = {0, 1, 4, 5};
    row.m_delta_weight = 1.0;
    row.m_dst_onv = {0, 1, 4, 5};
    row.m_delta_weight = 1.0;
}

TEST(Solver, Consolidation) {
    fciqmc_config::Document opts;
    opts.m_propagator.m_nw_target = 1000;
    opts.m_propagator.m_consolidate_spawns = true;
    opts.m_wavefunction.m_nw_init = 10;
    opts.verify();
    FermionHamiltonian ham(defs::assets_root + "/HF_RDMs/FCIDUMP", false);
    Wavefunction wf(opts, ham.nsite());
    props::Exact prop(ham, opts, wf.m_format);
    buffered::FrmOnv ref_onv(ham.nsite());
    ham.set_hf_onv(ref_onv, 0);

    auto ref_energy = ham.get_energy(ref_onv);
    auto ref_loc = wf.create_row(0, ref_onv, ref_energy, 1);
    wf.set_weight(opts.m_wavefunction.m_nw_init);

    Solver solver(opts, prop, wf, ref_loc);

    /*
     * spoof a recvd list
     */

    solver.loop_over_spawned();

}