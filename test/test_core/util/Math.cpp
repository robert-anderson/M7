//
// Created by rja on 14/06/22.
//

#include "test_core/defs.h"
#include "M7_lib/util/Math.h"

using namespace math;

TEST(UtilMath, CompileTimePow) {
    ASSERT_EQ(pow<3>(5), 5 * 5 * 5);
    ASSERT_EQ(pow<10>(2), 1024);
    ASSERT_EQ(pow<0>(10), 1ul);
    ASSERT_EQ(pow<1>(10), 10ul);
    ASSERT_EQ(pow<1>(0), 0ul);
}

TEST(UtilMath, GeoMean) {
    v_t<double> v = {0.12 , 0.009, 0.034, 0.074, 0.023, 0.003, 0.002, 0.079};
    ASSERT_NEAR_EQ(math::geo_mean(v), 0.02036831352);
}

TEST(UtilMath, Logarithm) {
    const double d1 = -8;
    const double d2 = 4;
    ASSERT_NEAR_EQ((math::logarithm(d1) + math::logarithm(d2)).exp(), d1 * d2);
    ASSERT_NEAR_EQ((math::logarithm(d1) - math::logarithm(d2)).exp(), d1 / d2);
}