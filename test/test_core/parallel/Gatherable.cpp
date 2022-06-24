//
// Created by Robert J. Anderson on 12/07/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/Gatherable.h"

TEST(Gatherable, Test){
    Gatherable<uint_t> gatherable;
    gatherable = mpi::irank();
    auto result = gatherable.mpi_gather();
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank){
        ASSERT_EQ(irank, result[irank]);
    }
}