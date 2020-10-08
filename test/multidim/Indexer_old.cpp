//
// Created by Robert John Anderson on 2020-02-22.
//

#include <gtest/gtest.h>
#include "src/core/multidim/Indexer_old.h"

TEST(Indexer, CorrectMapping) {
    Indexer_old<3> indexer(2, 4, 3);
    auto shape = std::array<size_t, 3>{2, 4, 3};
    ASSERT_EQ(indexer.shape(), shape);
    auto strides = std::array<size_t, 3>{12, 3, 1};
    ASSERT_EQ(indexer.strides(), strides);
    ASSERT_EQ(indexer.get(0, 0, 0), 0);
    ASSERT_EQ(indexer.get(0, 0, 1), 1);
    ASSERT_EQ(indexer.get(0, 0, 2), 2);
    ASSERT_EQ(indexer.get(0, 1, 0), 3);
    ASSERT_EQ(indexer.get(0, 1, 1), 4);
    ASSERT_EQ(indexer.get(0, 1, 2), 5);
    ASSERT_EQ(indexer.get(0, 2, 0), 6);
    ASSERT_EQ(indexer.get(0, 2, 1), 7);
    ASSERT_EQ(indexer.get(0, 2, 2), 8);
    ASSERT_EQ(indexer.get(0, 3, 0), 9);
    ASSERT_EQ(indexer.get(0, 3, 1), 10);
    ASSERT_EQ(indexer.get(0, 3, 2), 11);
    ASSERT_EQ(indexer.get(1, 0, 0), 12);
    ASSERT_EQ(indexer.get(1, 0, 1), 13);
    ASSERT_EQ(indexer.get(1, 0, 2), 14);
    ASSERT_EQ(indexer.get(1, 1, 0), 15);
    ASSERT_EQ(indexer.get(1, 1, 1), 16);
    ASSERT_EQ(indexer.get(1, 1, 2), 17);
    ASSERT_EQ(indexer.get(1, 2, 0), 18);
    ASSERT_EQ(indexer.get(1, 2, 1), 19);
    ASSERT_EQ(indexer.get(1, 2, 2), 20);
    ASSERT_EQ(indexer.get(1, 3, 0), 21);
    ASSERT_EQ(indexer.get(1, 3, 1), 22);
    ASSERT_EQ(indexer.get(1, 3, 2), 23);

    ASSERT_EQ(indexer.nelement(), 24);
}


TEST(Indexer, ReducedIndexSubscripting) {
    Indexer_old<4> indexer(5, 2, 4, 3);
    ASSERT_EQ(indexer.get_sub(), 0);
    ASSERT_EQ(indexer.get_sub(0), 0);
    ASSERT_EQ(indexer.get_sub(1), 1 * (2 * 4 * 3));
    ASSERT_EQ(indexer.get_sub(1, 2), 1 * (2 * 4 * 3) + 2 * (4 * 3));
    ASSERT_EQ(indexer.get_sub(1, 1, 1), 1 * (2 * 4 * 3) + 1 * (4 * 3) + 1 * (3));
}