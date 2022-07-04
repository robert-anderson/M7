//
// Created by Robert J. Anderson on 13/07/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/SharedArray.h"

TEST(SharedArray, Test) {
    const uint_t nelement_per_rank = 10000;
    const uint_t nelement = nelement_per_rank*mpi::nrank();
    SharedArray<double> array(nelement);

    if (mpi::on_node_i_am_root()) array.set(nelement-1, 1234);
    mpi::barrier_on_node();
    std::cout << array[nelement-1] <<std::endl;
    ASSERT(array[nelement-1] == 1234)
//    if (mpi::on_node_i_am_root()) {
//        for (uint_t i = 0ul; i < array.size(); ++i) {
//            array.set(i, i);
//        }
//    }
//    ASSERT(array[nelement_per_rank-1]==nelement_per_rank-1)
//    for (uint_t i = 0ul; i < array.size(); ++i) {
//        ASSERT_EQ(array[i], i);
//    }
}

TEST(SharedArray, VectorTest) {
    const uint_t nrow = 10;
    const uint_t nelement_per_rank = 10;
    v_t<SharedArray<double>> arrays{};
    for (uint_t i = 0ul; i < nrow; ++i) {
        arrays.emplace_back(mpi::nrank() * nelement_per_rank);
    }
}
