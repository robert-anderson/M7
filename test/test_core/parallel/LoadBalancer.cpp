//
// Created by rja on 18/07/22.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/LoadBalancer.h"
#include "M7_lib/util/Hash.h"

TEST(LoadBalancer, RankReallocator) {
    const uint_t nrank = 10;
    const uint_t nblock = 123;
    const auto block_iranks = convert::vector<uint_t>(hash::in_range(0, nblock, 0, nrank));
    const auto work_figs = convert::vector<double>(hash::in_range(1, nblock, 4, 23));

    RankReallocator rr(block_iranks, work_figs, nrank);
    const uintv_t chk_moving_blocks = {44, 76, 98, 44, 98, 76, 44, 18, 1, 53, 42, 36, 21};
    ASSERT_EQ(rr.m_moving_iblocks, chk_moving_blocks);
    const uintv_t chk_dst_iranks = {6, 6, 9, 3, 4, 6, 9, 3, 1, 4, 5, 6, 9};
    ASSERT_EQ(rr.m_dst_iranks, chk_dst_iranks);
}