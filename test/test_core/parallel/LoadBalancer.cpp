//
// Created by rja on 18/07/22.
//

#include "gtest/gtest.h"
#include "M7_lib/communication/Communicator.h"
#include "M7_lib/util/Hash.h"

TEST(LoadBalancer, Redistributor) {
    const uint_t nrank = 10;
    const uint_t nblock = 123;
    const auto block_iranks = convert::vector<uint_t>(hash::in_range(0, nblock, 0, nrank));
    const auto work_figs = convert::vector<double>(hash::in_range(1, nblock, 4, 23));

    Redistributor rr(block_iranks, work_figs, nrank);
    const v_t<Redistributor::Move> chk_moves = {
            {12, 4}, {52, 4}, {13, 4}, {19, 4}, {53, 4}, {66, 1}, {103, 3},
            {9, 4}, {33, 1}, {39, 5}, {57, 3}, {11, 4}, {21, 6}, {116, 4}
    };

    ASSERT_EQ(rr.m_moves, chk_moves);
}