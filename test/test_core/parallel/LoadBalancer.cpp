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
            {44, 6}, {76, 6}, {98, 9}, {44, 3}, {98, 4}, {76, 6},
            {44, 9}, {18, 3}, {1, 1}, {53, 4}, {42, 5}, {36, 6}, {21, 9}
    };
    ASSERT_EQ(rr.m_moves, chk_moves);
}