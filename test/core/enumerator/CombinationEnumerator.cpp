//
// Created by Robert John Anderson on 2020-03-30.
//

#include <src/core/enumerator/CombinationEnumerator.h>
#include "gtest/gtest.h"

std::vector<std::vector<size_t>> get_results() {
    return {
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
}

TEST(CombinationEnumerator, Iterate) {

    CombinationEnumerator enumerator(6, 2);
    defs::inds inds(2);
    for (auto result : get_results()) {
        enumerator.next(inds);
        ASSERT_EQ(inds[0], result[0]);
        ASSERT_EQ(inds[1], result[1]);
    }
}

TEST(CombinationEnumerator, NextMethod) {
    CombinationEnumerator enumerator(6, 2);
    defs::inds inds(2);
    auto results = get_results();
    size_t i = ~0ul;
    while (enumerator.next(inds, i)) {
        ASSERT_EQ(inds, results[i]);
    }
    ASSERT_EQ(i, results.size());
}