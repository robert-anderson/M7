//
// Created by rja on 18/05/2021.
//

#include <src/core/dynamics/Wavefunction.h>
#include <src/core/basis/CiSpaces.h>
#include "gtest/gtest.h"

TEST(Subspace, Test){
    fciqmc_config::Document opts;
    opts.m_propagator.m_nw_target = 1000;
    opts.m_wavefunction.m_load_balancing.m_nblock_per_rank = 3;
    opts.verify();
    const size_t nsite = 6;
    const size_t nelec = 6;
    Wavefunction wf(opts, nsite);
    ci_gen::SpinSym gen(nsite, nelec, 0, ci_gen::default_include_fn(wf));

    auto& table = wf.m_store;
    gen(table.m_row, table.m_row.m_onv);
    ASSERT_TRUE(wf.m_ra.verify());

    /*
     * select some arbitrary number of the rows from the beginning of the walker table
     */
    const size_t nrow_rank_lo = 4, nrow_rank_hi = 7;
    ASSERT_LT(nrow_rank_lo, table.m_hwm);
    ASSERT_LT(nrow_rank_hi, table.m_hwm);
    const size_t nrow_this_rank = hashing::in_range(mpi::irank(), nrow_rank_lo, nrow_rank_hi);


    Wavefunction::DynamicRowSet rowset(wf, "some row set");
    /*
     * and add them to the local subspace
     */
    auto row = table.m_row;
    for (row.restart(); row.in_range(nrow_this_rank); row.step()){
        rowset.add_(row);
    }
}