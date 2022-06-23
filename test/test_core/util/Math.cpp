//
// Created by rja on 14/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Math.h"

using namespace math;

TEST(UtilMath, CompileTimePow) {
    ASSERT_EQ(pow<3>(5), 5 * 5 * 5);
    ASSERT_EQ(pow<10>(2), 1024);
    ASSERT_EQ(pow<0>(10), 1ul);
    ASSERT_EQ(pow<1>(10), 10ul);
    ASSERT_EQ(pow<1>(0), 0ul);
}