//
// Created by rja on 06/02/2021.
//

#include "gtest/gtest.h"
#include "src/core/nd/NdSequence.h"

TEST(NdSequence, Test) {
    const size_t n0 = 3, n1 = 4, n2 = 6;
    NdSequence<3> sequence({n0, n1, n2});
    auto cursor = sequence.cursor();

    for (size_t i0 = 0ul; i0 < n0; ++i0) {
        for (size_t i1 = 0ul; i1 < n1; ++i1) {
            for (size_t i2 = 0ul; i2 < n2; ++i2) {
                ASSERT_EQ(i0, cursor.inds()[0]);
                ASSERT_EQ(i1, cursor.inds()[1]);
                ASSERT_EQ(i2, cursor.inds()[2]);
                cursor++;
            }
        }
    }
}