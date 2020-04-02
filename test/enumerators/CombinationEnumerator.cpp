//
// Created by Robert John Anderson on 2020-02-16.
//

#include <gtest/gtest.h>

#if 0
TEST(CombinationEnumerator, Test) {
    std::vector<std::vector<size_t>> results =
        {
            {0, 1},
            {0, 2},
            {0, 3},
            {0, 4},
            {0, 5},
            {1, 2},
            {1, 3},
            {1, 4},
            {1, 5},
            {2, 3},
            {2, 4},
            {2, 5},
            {3, 4},
            {3, 5},
            {4, 5}
        };
    CombinationEnumerator enumerator(6, 2);
    defs::inds inds(2);
    for (auto result : results) {
        enumerator.next(inds);
        ASSERT_EQ(inds[0], result[0]);
        ASSERT_EQ(inds[1], result[1]);
    }
}
#endif