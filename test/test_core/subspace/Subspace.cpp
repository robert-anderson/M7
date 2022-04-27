//
// Created by Robert J. Anderson on 18/05/2021.
//

#include <M7_lib/wavefunction/Wavefunction.h>
#include "gtest/gtest.h"

#if 0
TEST(Subspace, Test){
    conf::Document opts;
    opts.m_propagator.m_nw_target = 1000;
    opts.m_wavefunction.m_load_balancing.m_nblock_per_rank = 3;
    opts.verify();
    const size_t nsite = 6;
    const size_t nelec = 6;
    const BasisData bd = {nsite, {}};
    Wavefunction wf(opts, bd);
    ci_gen::SpinSym gen(bd, nelec, 0, ci_gen::default_include_fn(wf));

    auto& table = wf.m_store;
    gen(table.m_row, table.m_row.m_mbf);
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
#endif