//
// Created by Robert John Anderson on 2020-02-22.
//

#include <gtest/gtest.h>
#include "src/multidim/NdArray.h"

TEST(NdArray, CorrectMapping) {
    NdArray<double, 3> array(2, 4, 3);
    auto shape = std::array<size_t, 3>{2, 4, 3};
    ASSERT_EQ(array.shape(), shape);
    auto strides = std::array<size_t, 3>{12, 3, 1};
    ASSERT_EQ(array.strides(), strides);

    *array.view(1, 2, 2) = 4.5;
    ASSERT_EQ(*array.view(1, 2, 2), 4.5);
}

TEST(NdArray, SubArrays) {
    NdArray<double, 3> array(2, 4, 3);
    auto subarray = array.subarray(1);
}