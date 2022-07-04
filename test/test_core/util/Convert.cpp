//
// Created by rja on 16/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Convert.h"

TEST(UtilConvert, IntVectorToString){
    v_t<int> v = {1, 5, 3, 6};
    ASSERT_EQ(convert::to_string(v), "[1, 5, 3, 6]");
}

TEST(UtilConvert, IntVectorVectorToString){
    v_t<v_t<int>> v = {{1, 4}, {5}, {4, 5, 3}, {5, 6, 6, 0}};
    ASSERT_EQ(convert::to_string(v), "[[1, 4], [5], [4, 5, 3], [5, 6, 6, 0]]");
}

TEST(UtilConvert, FloatVectorToString) {
    v_t<float> v = {1.0/11.0, -1.234, 5.3, 3.01, 1.0/7.0};
    ASSERT_EQ(convert::to_string(v, 6), "[9.090909e-02, -1.234000e+00, 5.300000e+00, 3.010000e+00, 1.428571e-01]");
    ASSERT_EQ(convert::to_string(v, 1), "[9.1e-02, -1.2e+00, 5.3e+00, 3.0e+00, 1.4e-01]");
    ASSERT_EQ(convert::to_string(v, 0), "[9e-02, -1e+00, 5e+00, 3e+00, 1e-01]");
}