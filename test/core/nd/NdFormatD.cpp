//
// Created by rja on 30/03/2021.
//

#include "gtest/gtest.h"
#include "src/core/nd/NdFormatD.h"

TEST(NdFormatD, Test) {
    NdFormatD format({2, 4, 3});
    size_t i=~0ul;
    ASSERT_EQ(format.flatten({0,0,0}), ++i);
    ASSERT_EQ(format.flatten({0,0,1}), ++i);
    ASSERT_EQ(format.flatten({0,0,2}), ++i);
    ASSERT_EQ(format.flatten({0,1,0}), ++i);
    ASSERT_EQ(format.flatten({0,1,1}), ++i);
    ASSERT_EQ(format.flatten({0,1,2}), ++i);
    ASSERT_EQ(format.flatten({0,2,0}), ++i);
    ASSERT_EQ(format.flatten({0,2,1}), ++i);
    ASSERT_EQ(format.flatten({0,2,2}), ++i);
    ASSERT_EQ(format.flatten({0,3,0}), ++i);
    ASSERT_EQ(format.flatten({0,3,1}), ++i);
    ASSERT_EQ(format.flatten({0,3,2}), ++i);
    ASSERT_EQ(format.flatten({1,0,0}), ++i);
    ASSERT_EQ(format.flatten({1,0,1}), ++i);
    ASSERT_EQ(format.flatten({1,0,2}), ++i);
    ASSERT_EQ(format.flatten({1,1,0}), ++i);
    ASSERT_EQ(format.flatten({1,1,1}), ++i);
    ASSERT_EQ(format.flatten({1,1,2}), ++i);
    ASSERT_EQ(format.flatten({1,2,0}), ++i);
    ASSERT_EQ(format.flatten({1,2,1}), ++i);
    ASSERT_EQ(format.flatten({1,2,2}), ++i);
    ASSERT_EQ(format.flatten({1,3,0}), ++i);
    ASSERT_EQ(format.flatten({1,3,1}), ++i);
    ASSERT_EQ(format.flatten({1,3,2}), ++i);

}