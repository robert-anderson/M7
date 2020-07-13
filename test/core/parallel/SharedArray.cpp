//
// Created by rja on 13/07/2020.
//

#include "gtest/gtest.h"
#include "src/core/parallel/SharedArray.h"

TEST(SharedArray, Test){
    const size_t nelement_per_rank = 10;
    SharedArray<double> array(mpi::nrank()*nelement_per_rank);
    if (mpi::on_node_i_am_root()) {
        for (size_t i = 0ul; i < array.size(); ++i) {
            array[i] = i;
        }
    }
    for (size_t i = 0ul; i < array.size(); ++i) {
        ASSERT_EQ(array[i], i);
    }
}
