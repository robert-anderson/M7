//
// Created by rja on 15/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Array.h"

TEST(UtilArray, ToVector){
    std::array<float, 5> a = {4.5, -3.4, 3.5, 0.2, -5.6};
    auto v = array::to_vector(a);
    std::vector<float> v_chk = {4.5, -3.4, 3.5, 0.2, -5.6};
    ASSERT_EQ(v, v_chk);
}

TEST(UtilArray, Filled){
    auto a = array::filled<float, 5>(6.7);
    std::array<float, 5> a_chk = {6.7, 6.7, 6.7, 6.7, 6.7};
    ASSERT_EQ(a, a_chk);
}