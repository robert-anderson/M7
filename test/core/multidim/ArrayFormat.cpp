//
// Created by RJA on 18/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/multidim/ArrayFormat.h"

TEST(ArrayFormat, OneDimCase){
    ArrayFormat<1> af(12);
    for (size_t i = 0ul; i < af.nelement(); ++i) {
        ASSERT_EQ(af[i], i);
    }
}

TEST(ArrayFormat, Sequence) {
    ArrayFormat<3> af(4, 2, 3);
    size_t i = 0ul;
    for (size_t i0 = 0ul; i0 < af.shape()[0]; ++i0) {
        for (size_t i1 = 0ul; i1 < af.shape()[1]; ++i1) {
            for (size_t i2 = 0ul; i2 < af.shape()[2]; ++i2) {
                ASSERT_EQ(af[i0][i1][i2], i);
                ++i;
            }
        }
    }
}