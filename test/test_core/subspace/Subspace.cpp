//
// Created by Robert J. Anderson on 18/05/2021.
//

#include <M7_lib/wavefunction/Wavefunction.h>
#include "gtest/gtest.h"

#if 0
TEST(Subspace, Test){
    conf::Document opts;
    opts.m_propagator.m_nw_target = 1000;
    opts.m_wavefunction.m_distribution.m_nblock_per_rank = 3;
    opts.verify();
    const uint_t nsite = 6;
    const uint_t nelec = 6;
    const BasisData bd = {nsite, {}};
    Fci wf(opts, bd);
    ci_gen::SpinSym gen(bd, nelec, 0, ci_gen::default_include_fn(wf));

    auto& table = wf.m_store;
    gen(table.m_row, table.m_row.m_mbf);
    ASSERT_TRUE(wf.m_ra.verify());

    /*
     * select some arbitrary number of the rows from the beginning of the walker table
     */
    const uint_t nrow_rank_lo = 4, nrow_rank_hi = 7;
    ASSERT_LT(nrow_rank_lo, table.nrow_in_use());
    ASSERT_LT(nrow_rank_hi, table.nrow_in_use());
    const uint_t nrow_this_rank = hashing::in_range(mpi::irank(), nrow_rank_lo, nrow_rank_hi);


    Fci::DynamicRowSet rowset(wf, "some row set");
    /*
     * and add them to the local subspace
     */
    auto row = table.m_row;
    for (row.restart(); row.in_range(nrow_this_rank); row.step()){
        rowset.add_(row);
    }
}
#endif