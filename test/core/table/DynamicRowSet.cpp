//
// Created by rja on 19/01/2021.
//

#include "gtest/gtest.h"
#include "src/core/field/Elements.h"
#include "src/core/parallel/RankAllocator.h"
#include "TableTest.h"


TEST(DynamicRowSet, Test) {
    const size_t nsite = 12;
    const size_t nrow_per_rank = 4;
    BufferedTable<table_test::DetMappedTable> bt("Test table", nsite);
    RankAllocator<fields::Det> ra(100, 10);

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

    for (size_t ibit = 2; ibit < 2 * nsite - 2; ++ibit) {
        det.set(ibit);
        if (bt.m_hwm < nrow_per_rank and mpi::i_am(ra.get_rank(det))) {
            auto irow = bt.push_back();
            bt.m_config(irow) = det;
        }
        det.clr(ibit);
    }
    bt.print_contents();

    auto drs = bt.dynamic_row_set(ra);

    drs.add(0);
}