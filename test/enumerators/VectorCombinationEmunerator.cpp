//
// Created by Robert John Anderson on 2020-02-16.
//

#if 0
#include "gtest/gtest.h"
#include <src/enumerators/VectorCombinationEnumerator.h>
#include <src/utils.h>

TEST(VectorCombinationEnumerator, Test) {
    std::vector<std::vector<size_t>> results =
        {
            {10, 11},
            {10, 12},
            {10, 13},
            {10, 14},
            {10, 15},
            {11, 12},
            {11, 13},
            {11, 14},
            {11, 15},
            {12, 13},
            {12, 14},
            {12, 15},
            {13, 14},
            {13, 15},
            {14, 15}
        };
    defs::inds vec{10, 11, 12, 13, 14, 15};
    VectorCombinationEnumerator enumerator(vec, 2);
    defs::inds inds(2);
    for (auto result : results) {
        enumerator.next(inds);
        ASSERT_EQ(inds[0], result[0]);
        ASSERT_EQ(inds[1], result[1]);
    }
}

TEST(VectorCombinationEnumerator, Nested) {
    defs::inds outer{10, 11, 12, 13, 14};
    defs::inds inner{20, 22, 24};
    VectorCombinationEnumerator outer_enumerator(outer, 3);
    defs::inds outer_inds(3);
    while(outer_enumerator.next(outer_inds)){
        {
            VectorCombinationEnumerator inner_enumerator(inner, 2);
            defs::inds inner_inds(2);
            while(inner_enumerator.next(inner_inds)) {
                std::cout << utils::to_string(outer_inds) << utils::to_string(inner_inds) << std::endl;
            }
        }
    }
}
#endif