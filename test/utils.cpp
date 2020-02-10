//
// Created by Robert John Anderson on 2020-01-18.
//

#include <gtest/gtest.h>
#include "../src/utils.h"

using namespace integer_utils;

TEST(utils, factorial){
    ASSERT_DEATH(factorial(-1), R"(Assertion.*)");
    ASSERT_EQ(factorial(0), 1ul);
    ASSERT_EQ(factorial(1), 1ul);
    ASSERT_EQ(factorial(2), 2ul);
    ASSERT_EQ(factorial(3), 6ul);
    ASSERT_EQ(factorial(20), 2432902008176640000ul);
    // overflows at 21!
    ASSERT_EQ(factorial(21), 14197454024290336768ul);
}

TEST(utils, combinatorial){
    for (auto n{0ul}; n<20; ++n) {
        for (auto r{0ul}; r <= n; ++r) {
            ASSERT_EQ(combinatorial(n, r),
                      factorial(n) / (factorial(n - r) * factorial(r)));
        }
    }
    ASSERT_EQ(combinatorial(100, 5), 75287520ul);
    ASSERT_EQ(combinatorial(100, 10), 17310309456440ul);
}
