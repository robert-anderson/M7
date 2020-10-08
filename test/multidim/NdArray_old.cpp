//
// Created by Robert John Anderson on 2020-02-22.
//

#include <gtest/gtest.h>
#include "src/core/multidim/NdArray.h"

TEST(NdArray, CorrectMapping) {
//    NdArray<double, 3> array(2, 4, 3);
//    auto shape = std::array<size_t, 3>{2, 4, 3};
//    ASSERT_EQ(array.shape(), shape);
//    auto strides = std::array<size_t, 3>{12, 3, 1};
//    ASSERT_EQ(array.strides(), strides);
//
//    *array.view(1, 2, 2) = 4.5;
//    ASSERT_EQ(*array.view(1, 2, 2), 4.5);
}

/*
TEST(NdArray, SubArrays) {
    NdArray<double, 4> array(2, 4, 3, 7);
    auto subarray = array.subarray(1, 2);
    auto subarray_shape = std::array<size_t, 2>{3, 7};
    ASSERT_EQ(subarray_shape, subarray.shape());
    *subarray.view(4, 3) = 1234.567;
    ASSERT_EQ(*array.view(1, 2, 4, 3), 1234.567);
}*/

TEST(NdArray, VectorCase) {
//    NdArray<double, 1> array(8);
//    std::vector<double> v(8, 8.9);
//    for (size_t i=0ul; i<v.size(); ++i) v[i]*=i;
//    array = v;
//    for (size_t i=0ul; i<v.size(); ++i) {
//        ASSERT_EQ(*array.view(i), v[i]);
//    }
}