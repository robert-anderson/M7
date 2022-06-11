//
// Created by rja on 04/06/22.
//

#include "test_core/defs.h"
#include "M7_lib/integrals/IntegralArray1e.h"


TEST(IntegralArray1e, SymNone_real) {
    typedef double T;
    integrals_1e::SymNone<T> ints(6);
    ASSERT_TRUE(ints.set(2, 5, 0.123));
    ASSERT_TRUE(ints.set(2, 5, 0.123));
    ASSERT_FALSE(ints.set(2, 5, 0.1234));
    ASSERT_COMPLEX_EQ(ints.get(2, 5), 0.123);
    ASSERT_COMPLEX_EQ(ints.get(5, 2), 0.0);
    ASSERT_TRUE(ints.set(5, 2, 0.234));
    ASSERT_COMPLEX_EQ(ints.get(5, 2), 0.234);
}

TEST(IntegralArray1e, SymNone_complex) {
    typedef std::complex<double> T;
    integrals_1e::SymNone<T> ints(6);
    const T v = {0.123, -0.234};
    ASSERT_TRUE(ints.set(2, 5, v));
    ASSERT_TRUE(ints.set(2, 5, v));
    ASSERT_FALSE(ints.set(2, 5, {0.1234, 0.2345}));
    ASSERT_COMPLEX_EQ(ints.get(2, 5), v);
    ASSERT_COMPLEX_EQ(ints.get(5, 2), 0.0);
}

TEST(IntegralArray1e, SymH_real) {
    typedef double T;
    integrals_1e::SymH<T> ints(6);
    ASSERT_TRUE(ints.set(2, 5, 0.123));
    ASSERT_TRUE(ints.set(2, 5, 0.123));
    ASSERT_FALSE(ints.set(2, 5, 0.1234));
    ASSERT_COMPLEX_EQ(ints.get(2, 5), 0.123);
    ASSERT_COMPLEX_EQ(ints.get(5, 2), 0.123);
    ASSERT_FALSE(ints.set(5, 2, 0.234));
    ASSERT_COMPLEX_EQ(ints.get(5, 2), 0.123);
}

TEST(IntegralArray1e, SymH_complex) {
    typedef std::complex<double> T;
    integrals_1e::SymH<T> ints(6);
    const T v = {0.123, -0.234};
    const T vc = {0.123, 0.234};
    ASSERT_TRUE(ints.set(2, 5, v));
    ASSERT_TRUE(ints.set(2, 5, v));
    ASSERT_FALSE(ints.set(2, 5, {0.1234, 0.2345}));
    ASSERT_COMPLEX_EQ(ints.get(2, 5), v);
    ASSERT_COMPLEX_EQ(ints.get(5, 2), vc);
    ASSERT_FALSE(ints.set(5, 2, v));
    ASSERT_TRUE(ints.set(5, 2, vc));
}