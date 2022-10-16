//
// Created by Robert J. Anderson on 19/01/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/RankAllocator.h"
#include "TableTest.h"

#if 0
TEST(DynamicRowSet, Test) {
    const uint_t nsite = 12;
    const uint_t nrow_per_rank = 4;
    BufferedTable<table_test::DetMappedTable> bt("Test table", nsite);
    RankAllocator<fields::Det> ra(bt, bt.m_key_field, 100, 10, 0.05);

    /*
     * arbitrary choice for testing purposes, but let's make the "determinants"
     * have the first two and last two bits set, and then we'll set another bit
     * until we have the required number of rows generated per rank.
     */
    elements::Det det(nsite);
    det.set(0);
    det.set(1);
    det.set(2 * nsite - 2);
    det.set(2 * nsite - 1);

    for (uint_t ibit = 2; ibit < 2 * nsite - 2; ++ibit) {
        det.set(ibit);
        if (bt.nrow_in_use() < nrow_per_rank and mpi::i_am(ra.get_rank(det))) {
            auto irow = bt.push_back();
            bt.m_config(irow) = det;
        }
        det.clr(ibit);
    }
    bt.print_contents();
}
#endif