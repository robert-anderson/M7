//
// Created by anderson on 03/06/2022.
//

#include "test_core/defs.h"
#include "M7_lib/integrals/IntegralArray2e.h"

TEST(IntegralArray2e_i, Sym2Fold_i_real) {
    typedef double T;
    IntegralArray2e_i<T> array(4);
    const std::array<T, 4> values = {0.123, 0.234, 0.345, 0.456};
    // -
    array.set(0, 1, 2, 3, values[0]);
    // r
    array.set(0, 1, 3, 2, values[1]);
    // rh
    array.set(2, 3, 0, 1, values[2]);
    // h
    array.set(3, 2, 0, 1, values[3]);

    // 0
    ASSERT_NEARLY_EQ(array.get(0, 1, 2, 3), values[0]);
    // i 0
    ASSERT_NEARLY_EQ(array.get(1, 0, 3, 2), values[0]);
    // 1
    ASSERT_NEARLY_EQ(array.get(0, 1, 3, 2), values[1]);
    // i 1
    ASSERT_NEARLY_EQ(array.get(1, 0, 2, 3), values[1]);
    // 2
    ASSERT_NEARLY_EQ(array.get(2, 3, 0, 1), values[2]);
    // i 2
    ASSERT_NEARLY_EQ(array.get(3, 2, 1, 0), values[2]);
    // 3
    ASSERT_NEARLY_EQ(array.get(3, 2, 0, 1), values[3]);
    // i 3
    ASSERT_NEARLY_EQ(array.get(2, 3, 1, 0), values[3]);
}

TEST(IntegralArray2e_i, Sym2Fold_i_complex) {
    typedef std::complex<double> T;
    IntegralArray2e_i<T> array(4);
    const std::array<T, 4> values = {T{0.123, 9.0}, T{0.234, -2.0}, T{0.345, -1.0}, T{0.456, 7.0}};
    // -
    array.set(0, 1, 2, 3, values[0]);
    // r
    array.set(0, 1, 3, 2, values[1]);
    // rh
    array.set(2, 3, 0, 1, values[2]);
    // h
    array.set(3, 2, 0, 1, values[3]);

    // 0
    ASSERT_NEARLY_EQ(array.get(0, 1, 2, 3), values[0]);
    // i 0
    ASSERT_NEARLY_EQ(array.get(1, 0, 3, 2), values[0]);
    // 1
    ASSERT_NEARLY_EQ(array.get(0, 1, 3, 2), values[1]);
    // i 1
    ASSERT_NEARLY_EQ(array.get(1, 0, 2, 3), values[1]);
    // 2
    ASSERT_NEARLY_EQ(array.get(2, 3, 0, 1), values[2]);
    // i 2
    ASSERT_NEARLY_EQ(array.get(3, 2, 1, 0), values[2]);
    // 3
    ASSERT_NEARLY_EQ(array.get(3, 2, 0, 1), values[3]);
    // i 3
    ASSERT_NEARLY_EQ(array.get(2, 3, 1, 0), values[3]);
}

TEST(IntegralArray2e_ih, Sym4Fold_ih_real) {
    typedef double T;
    IntegralArray2e_i<T> array(4);
    const std::array<T, 2> values = {0.123, 0.234};
    // -
    array.set(0, 1, 2, 3, values[0]);
    // r
    array.set(0, 1, 3, 2, values[1]);

    // 0
    ASSERT_NEARLY_EQ(array.get(0, 1, 2, 3), values[0]);
    // i 0
    ASSERT_NEARLY_EQ(array.get(1, 0, 3, 2), values[0]);
    // 1
    ASSERT_NEARLY_EQ(array.get(0, 1, 3, 2), values[1]);
    // i 1
    ASSERT_NEARLY_EQ(array.get(1, 0, 2, 3), values[1]);
    // h 0
    ASSERT_NEARLY_EQ(array.get(2, 3, 0, 1), values[0]);
    // ih 0
    ASSERT_NEARLY_EQ(array.get(3, 2, 1, 0), values[0]);
    // h 1
    ASSERT_NEARLY_EQ(array.get(3, 2, 0, 1), values[1]);
    // ih 1
    ASSERT_NEARLY_EQ(array.get(2, 3, 1, 0), values[1]);
}



TEST(IntegralArray2e_ihr, Sym8Fold_real) {
    IntegralArray2e_ihr<double> array(4);
    const std::array<defs::ham_t, 1> values = {0.123};
    // -
    array.set(0, 1, 2, 3, values[0]);

    ASSERT_NEARLY_EQ(array.get(0, 1, 2, 3), values[0]);
    ASSERT_NEARLY_EQ(array.get(1, 0, 3, 2), values[0]);
    ASSERT_NEARLY_EQ(array.get(0, 1, 3, 2), values[0]);
    ASSERT_NEARLY_EQ(array.get(1, 0, 2, 3), values[0]);
    ASSERT_NEARLY_EQ(array.get(2, 3, 0, 1), values[0]);
    ASSERT_NEARLY_EQ(array.get(3, 2, 1, 0), values[0]);
    ASSERT_NEARLY_EQ(array.get(3, 2, 0, 1), values[0]);
    ASSERT_NEARLY_EQ(array.get(2, 3, 1, 0), values[0]);
}