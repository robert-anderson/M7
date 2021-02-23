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
    opts.nwalker_initial = 10.0;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.01;
    opts.load_balance_period = 5;
    opts.nload_balance_block_per_rank = 40;
    opts.nwalker_target = 100000;
    opts.shift_damp = 0.5;
    opts.shift_initial = 0.0;
    opts.ncycle = 500000;
    opts.init();
//const auto benchmark = -108.81138657563143;
    FermionHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv ref_onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        ref_onv.set(0, i);
        ref_onv.set(1, i);
    }
    StochasticPropagator prop(ham, opts);
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    ASSERT_EQ(&wf.m_store.m_row.m_onv, &KeyField<WalkerTableRow>::get(wf.m_store.m_row));
    //const size_t store_nrow = (opts.walker_buffer_size_factor_initial*opts.nwalker_target)/mpi::nrank();
    //ASSERT_EQ(wf.m_store.m_nrow, store_nrow);
    //const size_t send_nrow = (opts.spawn_buffer_size_factor_initial*opts.nwalker_target)/mpi::nrank();
    //ASSERT_EQ(wf.m_comm.send().nrow_per_table(), send_nrow);

    auto ref_energy = ham.get_energy(ref_onv);

    auto ref_loc = wf.create_walker(ref_onv, opts.nwalker_initial, ref_energy, 1);
    prop.m_shift = ref_energy+opts.shift_initial;
    Solver solver(prop, wf, ref_loc);
    solver.execute(opts.ncycle);
}


#else

TEST(StochasticPropagator, BosonTest) {
    Options opts;
    opts.nwalker_initial = 10;
    opts.nadd_initiator = 3.0;
    opts.tau_initial = 0.01;
    opts.nwalker_target = 10000;
    opts.shift_damp = 0.4;
    opts.ncycle = 3000;

    // -10.328242246088791
    Hamiltonian<> ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 2, 1.4, 0.3);

    ASSERT_TRUE(ham.spin_conserving());
    buffered::Onv<> onv(ham.nsite());
    for (size_t i = 0ul; i < ham.nelec() / 2; ++i) {
        onv.m_fonv.set(0, i);
        onv.m_fonv.set(1, i);
    }
    Wavefunction wf(opts, ham.nsite());
    wf.m_store.expand(10);
    wf.m_comm.expand(800);
    StochasticPropagator prop(ham, opts);
    auto ref_energy = ham.get_energy(onv);
    prop.m_shift = ref_energy;//benchmark;
    auto ref_loc = wf.create_walker(onv, opts.nwalker_initial, ref_energy, 1);
    Solver solver(prop, wf, ref_loc);

    std::cout << "Reference Energy: " << ref_energy << std::endl;

    for (size_t i = 0ul; i < opts.ncycle; ++i) {
        solver.execute();
        std::cout << i << " " << wf.m_store.m_hwm << " " << std::sqrt(wf.square_norm()) << std::endl;
    }
}
#endif