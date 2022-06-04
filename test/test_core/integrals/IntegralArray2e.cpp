//
// Created by anderson on 03/06/2022.
//

#include "test_core/defs.h"
#include "M7_lib/integrals/IntegralArray2e.h"


TEST(IntegralArray2e, SymNone_real) {
    typedef double T;
    integrals_2e::SymNone<T> ints(6);
    const std::array<T, 8> v = {0.123, 0.234, 0.345, 0.456, 0.567, 0.678, 0.789, 0.890};
    // -
    ASSERT_TRUE(ints.set(0, 1, 2, 3, v[0]));
    // I
    ASSERT_TRUE(ints.set(1, 0, 3, 2, v[1]));
    // R
    ASSERT_TRUE(ints.set(0, 1, 3, 2, v[2]));
    // RI
    ASSERT_TRUE(ints.set(1, 0, 2, 3, v[3]));
    // H
    ASSERT_TRUE(ints.set(2, 3, 0, 1, v[4]));
    // HI
    ASSERT_TRUE(ints.set(3, 2, 1, 0, v[5]));
    // HR
    ASSERT_TRUE(ints.set(3, 2, 0, 1, v[6]));
    // HRI
    ASSERT_TRUE(ints.set(2, 3, 1, 0, v[7]));

    // -
    ASSERT_NEARLY_EQ(ints.get(0, 1, 2, 3), v[0]);
    // I
    ASSERT_NEARLY_EQ(ints.get(1, 0, 3, 2), v[1]);
    // R
    ASSERT_NEARLY_EQ(ints.get(0, 1, 3, 2), v[2]);
    // RI
    ASSERT_NEARLY_EQ(ints.get(1, 0, 2, 3), v[3]);
    // H
    ASSERT_NEARLY_EQ(ints.get(2, 3, 0, 1), v[4]);
    // HI
    ASSERT_NEARLY_EQ(ints.get(3, 2, 1, 0), v[5]);
    // HR
    ASSERT_NEARLY_EQ(ints.get(3, 2, 0, 1), v[6]);
    // HRI
    ASSERT_NEARLY_EQ(ints.get(2, 3, 1, 0), v[7]);
}

TEST(IntegralArray2e, SymH_real) {
    typedef double T;
    integrals_2e::SymH<T> ints(6);
    const std::array<T, 4> v = {0.123, 0.234, 0.345, 0.456};
    // -
    ASSERT_TRUE(ints.set(0, 1, 2, 3, v[0]));
    // I
    ASSERT_TRUE(ints.set(1, 0, 3, 2, v[1]));
    // R
    ASSERT_TRUE(ints.set(0, 1, 3, 2, v[2]));
    // RI
    ASSERT_TRUE(ints.set(1, 0, 2, 3, v[3]));

    // -
    ASSERT_NEARLY_EQ(ints.get(0, 1, 2, 3), v[0]);
    // I
    ASSERT_NEARLY_EQ(ints.get(1, 0, 3, 2), v[1]);
    // R
    ASSERT_NEARLY_EQ(ints.get(0, 1, 3, 2), v[2]);
    // RI
    ASSERT_NEARLY_EQ(ints.get(1, 0, 2, 3), v[3]);
    // H
    ASSERT_NEARLY_EQ(ints.get(2, 3, 0, 1), v[0]);
    // HI
    ASSERT_NEARLY_EQ(ints.get(3, 2, 1, 0), v[1]);
    // HR
    ASSERT_NEARLY_EQ(ints.get(3, 2, 0, 1), v[2]);
    // HRI
    ASSERT_NEARLY_EQ(ints.get(2, 3, 1, 0), v[3]);
}

TEST(IntegralArray2e, SymI_real) {
    typedef double T;
    integrals_2e::SymI<T> ints(6);
    const std::array<T, 4> v = {0.123, 0.234, 0.345, 0.456};
    // -
    ASSERT_TRUE(ints.set(0, 1, 2, 3, v[0]));
    // R
    ASSERT_TRUE(ints.set(0, 1, 3, 2, v[1]));
    // H
    ASSERT_TRUE(ints.set(2, 3, 0, 1, v[2]));
    // HR
    ASSERT_TRUE(ints.set(3, 2, 0, 1, v[3]));

    // -
    ASSERT_NEARLY_EQ(ints.get(0, 1, 2, 3), v[0]);
    // I
    ASSERT_NEARLY_EQ(ints.get(1, 0, 3, 2), v[0]);
    // R
    ASSERT_NEARLY_EQ(ints.get(0, 1, 3, 2), v[1]);
    // RI
    ASSERT_NEARLY_EQ(ints.get(1, 0, 2, 3), v[1]);
    // H
    ASSERT_NEARLY_EQ(ints.get(2, 3, 0, 1), v[2]);
    // HI
    ASSERT_NEARLY_EQ(ints.get(3, 2, 1, 0), v[2]);
    // HR
    ASSERT_NEARLY_EQ(ints.get(3, 2, 0, 1), v[3]);
    // HRI
    ASSERT_NEARLY_EQ(ints.get(2, 3, 1, 0), v[3]);
}

TEST(IntegralArray2e, SymDH_real) {
    typedef double T;
    integrals_2e::SymDH<T> ints(6);
    const std::array<T, 2> v = {0.123, 0.234};
    // -
    ASSERT_TRUE(ints.set(0, 1, 2, 3, v[0]));
    // R
    ASSERT_TRUE(ints.set(0, 1, 3, 2, v[1]));

    // -
    ASSERT_NEARLY_EQ(ints.get(0, 1, 2, 3), v[0]);
    // I
    ASSERT_NEARLY_EQ(ints.get(1, 0, 3, 2), v[0]);
    // R
    ASSERT_NEARLY_EQ(ints.get(0, 1, 3, 2), v[1]);
    // RI
    ASSERT_NEARLY_EQ(ints.get(1, 0, 2, 3), v[1]);
    // H
    ASSERT_NEARLY_EQ(ints.get(2, 3, 0, 1), v[0]);
    // HI
    ASSERT_NEARLY_EQ(ints.get(3, 2, 1, 0), v[0]);
    // HR
    ASSERT_NEARLY_EQ(ints.get(3, 2, 0, 1), v[1]);
    // HRI
    ASSERT_NEARLY_EQ(ints.get(2, 3, 1, 0), v[1]);
}

TEST(IntegralArray2e, SymDR_real) {
    typedef double T;
    integrals_2e::SymDR<T> ints(6);
    const std::array<T, 2> v = {0.123, 0.234};
    // -
    ASSERT_TRUE(ints.set(0, 1, 2, 3, v[0]));
    // H
    ASSERT_TRUE(ints.set(2, 3, 0, 1, v[1]));

    // -
    ASSERT_NEARLY_EQ(ints.get(0, 1, 2, 3), v[0]);
    // I
    ASSERT_NEARLY_EQ(ints.get(1, 0, 3, 2), v[0]);
    // R
    ASSERT_NEARLY_EQ(ints.get(0, 1, 3, 2), v[0]);
    // RI
    ASSERT_NEARLY_EQ(ints.get(1, 0, 2, 3), v[0]);
    // H
    ASSERT_NEARLY_EQ(ints.get(2, 3, 0, 1), v[1]);
    // HI
    ASSERT_NEARLY_EQ(ints.get(3, 2, 1, 0), v[1]);
    // HR
    ASSERT_NEARLY_EQ(ints.get(3, 2, 0, 1), v[1]);
    // HRI
    ASSERT_NEARLY_EQ(ints.get(2, 3, 1, 0), v[1]);
}







TEST(IntegralArray2e, SymH_complex) {
    typedef std::complex<double> T;
    integrals_2e::SymH<T> ints(6);
    const std::array<T, 4> v = {T(0.123, -1), T(0.234, 2), T(0.345, -3), T(0.456, 4)};
    const std::array<T, 4> vc = {T(0.123, 1), T(0.234, -2), T(0.345, 3), T(0.456, -4)};
    // -
    ASSERT_TRUE(ints.set(0, 1, 2, 3, v[0]));
    // I
    ASSERT_TRUE(ints.set(1, 0, 3, 2, v[1]));
    // R
    ASSERT_TRUE(ints.set(0, 1, 3, 2, v[2]));
    // RI
    ASSERT_TRUE(ints.set(1, 0, 2, 3, v[3]));

    // -
    ASSERT_NEARLY_EQ(ints.get(0, 1, 2, 3), v[0]);
    // I
    ASSERT_NEARLY_EQ(ints.get(1, 0, 3, 2), v[1]);
    // R
    ASSERT_NEARLY_EQ(ints.get(0, 1, 3, 2), v[2]);
    // RI
    ASSERT_NEARLY_EQ(ints.get(1, 0, 2, 3), v[3]);
    // H
    ASSERT_NEARLY_EQ(ints.get(2, 3, 0, 1), vc[0]);
    // HI
    ASSERT_NEARLY_EQ(ints.get(3, 2, 1, 0), vc[1]);
    // HR
    ASSERT_NEARLY_EQ(ints.get(3, 2, 0, 1), vc[2]);
    // HRI
    ASSERT_NEARLY_EQ(ints.get(2, 3, 1, 0), vc[3]);
}