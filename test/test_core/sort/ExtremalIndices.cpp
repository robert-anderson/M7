//
// Created by Robert J. Anderson on 09/07/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/sort/ExtremalIndices.h"

TEST(ExtremalIndices, Ascending) {
    std::vector<int> data = {6, 13, -12, 2, 9, 0, 1, -99, 1999};
    auto cmp = [&data](const uint_t& i1, const uint_t& i2){
        return data[i1]<data[i2];
    };
    ExtremalIndices xi(cmp);
    xi.reset(data.size());
    xi.find(5);
    ASSERT_EQ(data[xi[0]], -99);
    ASSERT_EQ(data[xi[1]], -12);
    ASSERT_EQ(data[xi[2]], 0);
    ASSERT_EQ(data[xi[3]], 1);
    ASSERT_EQ(data[xi[4]], 2);
}

TEST(ExtremalIndices, AscendingAbs) {
    std::vector<int> data = {6, 13, -12, 2, 9, 0, 1, -99, 1999};
    auto cmp = [&data](const uint_t& i1, const uint_t& i2){
        return std::abs(data[i1])<std::abs(data[i2]);
    };
    ExtremalIndices xi(cmp);
    xi.reset(data.size());
    xi.find(5);
    ASSERT_EQ(data[xi[0]], 0);
    ASSERT_EQ(data[xi[1]], 1);
    ASSERT_EQ(data[xi[2]], 2);
    ASSERT_EQ(data[xi[3]], 6);
    ASSERT_EQ(data[xi[4]], 9);
}

TEST(ExtremalIndices, Descending) {
    std::vector<int> data = {6, 13, -12, 2, 9, 0, 1, -99, 1999};
    auto cmp = [&data](const uint_t& i1, const uint_t& i2){
        return data[i1]>data[i2];
    };
    ExtremalIndices xi(cmp);
    xi.reset(data.size());
    xi.find(5);
    ASSERT_EQ(data[xi[0]], 1999);
    ASSERT_EQ(data[xi[1]], 13);
    ASSERT_EQ(data[xi[2]], 9);
    ASSERT_EQ(data[xi[3]], 6);
    ASSERT_EQ(data[xi[4]], 2);
}

TEST(ExtremalIndices, DescendingAbs) {
    std::vector<int> data = {6, 13, -12, 2, 9, 0, 1, -99, 1999};
    auto cmp = [&data](const uint_t& i1, const uint_t& i2){
        return std::abs(data[i1])>std::abs(data[i2]);
    };
    ExtremalIndices xi(cmp);
    xi.reset(data.size());
    xi.find(5);
    ASSERT_EQ(data[xi[0]], 1999);
    ASSERT_EQ(data[xi[1]], -99);
    ASSERT_EQ(data[xi[2]], 13);
    ASSERT_EQ(data[xi[3]], -12);
    ASSERT_EQ(data[xi[4]], 9);
}

TEST(ExtremalIndices, AscendingMultiple) {
    // find elements in three separate calls, growing the sorted set each time
    std::vector<int> data = {6, 13, -12, 2, 9, 0, 1, -99, 1999};
    auto cmp = [&data](const uint_t& i1, const uint_t& i2){
        return data[i1]<data[i2];
    };
    ExtremalIndices xi(cmp);
    xi.reset(data.size());
    xi.find(3);
    ASSERT_EQ(xi.nfound(), 3);
    ASSERT_EQ(data[xi[0]], -99);
    ASSERT_EQ(data[xi[1]], -12);
    ASSERT_EQ(data[xi[2]], 0);
    xi.find(2);
    ASSERT_EQ(xi.nfound(), 5);
    ASSERT_EQ(data[xi[0]], -99);
    ASSERT_EQ(data[xi[1]], -12);
    ASSERT_EQ(data[xi[2]], 0);
    ASSERT_EQ(data[xi[3]], 1);
    ASSERT_EQ(data[xi[4]], 2);
    xi.find(2);
    ASSERT_EQ(xi.nfound(), 7);
    ASSERT_EQ(data[xi[0]], -99);
    ASSERT_EQ(data[xi[1]], -12);
    ASSERT_EQ(data[xi[2]], 0);
    ASSERT_EQ(data[xi[3]], 1);
    ASSERT_EQ(data[xi[4]], 2);
    ASSERT_EQ(data[xi[5]], 6);
    ASSERT_EQ(data[xi[6]], 9);
}