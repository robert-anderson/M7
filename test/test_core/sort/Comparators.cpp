//
// Created by rja on 13/07/2021.
//

#include "M7_lib/sort/Comparators.h"
#include "gtest/gtest.h"

TEST(Comparators, ValueFnsInt) {
    using namespace comparators;
    auto cmp_fn = make_value_cmp_fn<int>(false, false);
    // 1 is a better value than 2 wrt this comparator
    ASSERT_TRUE(cmp_fn(1, 2));
    // all comparator fns emitted by make_value_cmp_fn are strict, so 1 is not a better value than itself
    ASSERT_FALSE(cmp_fn(1, 1));
    // not taking the absolute value so -1 is better than 1
    ASSERT_TRUE(cmp_fn(-1, 1));

    cmp_fn = make_value_cmp_fn<int>(true, false);
    // 1 is a better value than 2 wrt this comparator
    ASSERT_TRUE(cmp_fn(1, 2));
    ASSERT_FALSE(cmp_fn(1, 1));
    // taking the absolute value so -1 is better not better than 1
    ASSERT_FALSE(cmp_fn(-1, 1));
    // but -2 is
    ASSERT_FALSE(cmp_fn(-2, 1));

    // below are the same tests with opposite sense
    cmp_fn = make_value_cmp_fn<int>(false, true);
    ASSERT_FALSE(cmp_fn(1, 2));
    ASSERT_FALSE(cmp_fn(1, 1));
    ASSERT_FALSE(cmp_fn(-1, 1));

    cmp_fn = make_value_cmp_fn<int>(true, true);
    ASSERT_FALSE(cmp_fn(1, 2));
    ASSERT_FALSE(cmp_fn(1, 1));
    ASSERT_FALSE(cmp_fn(-1, 1));
    ASSERT_TRUE(cmp_fn(-2, 1));
}