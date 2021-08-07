//
// Created by rja on 07/08/2021.
//

#include "src/core/util/utils.h"
#include "gtest/gtest.h"

TEST(Utils, CompileTimePow){
    ASSERT_EQ(utils::pow<3>(5), 5*5*5);
    ASSERT_EQ(utils::pow<10>(2), 1024);
    ASSERT_EQ(utils::pow<0>(10), 1ul);
    ASSERT_EQ(utils::pow<1>(10), 10ul);
    ASSERT_EQ(utils::pow<1>(0), 0ul);
}
TEST(Utils, CompileTimeNtup){
    ASSERT_EQ(utils::ntup<4>(15), integer_utils::combinatorial(15, 4));
    ASSERT_EQ(utils::ntup<1>(15), 15);
    ASSERT_EQ(utils::ntup<1>(1), 1);
    ASSERT_EQ(utils::ntup<0>(1), 1);
    ASSERT_EQ(utils::ntup<0>(15), 1);
}