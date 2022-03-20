//
// Created by rja on 12/07/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/Gatherable.h"

TEST(Gatherable, Test){
    Gatherable<size_t> gatherable;
    gatherable = mpi::irank();
    auto result = gatherable.mpi_gather();
    for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
        ASSERT_EQ(irank, result[irank]);
    }
}