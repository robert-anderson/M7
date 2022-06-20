//
// Created by rja on 16/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Convert.h"

TEST(UtilConvert, IntVectorToString){
    std::vector<int> v = {1, 5, 3, 6};
    ASSERT_EQ(utils::convert::to_string(v), "[1, 5, 3, 6]");
}

TEST(UtilConvert, IntVectorVectorToString){
    std::vector<std::vector<int>> v = {{1, 4}, {5}, {4, 5, 3}, {5, 6, 6, 0}};
    ASSERT_EQ(utils::convert::to_string(v), "[[1, 4], [5], [4, 5, 3], [5, 6, 6, 0]]");
}