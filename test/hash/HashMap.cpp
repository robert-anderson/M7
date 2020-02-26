//
// Created by Robert John Anderson on 2020-02-20.
//


#include <gtest/gtest.h>
#include "src/hash/HashMap.h"


#include <gtest/gtest.h>
#include <src/hash/HashMap.h>

#if 0

TEST(HashMap, SimpleLookup) {
    const size_t n = 100;
    HashMap <size_t> map(n);
    map.insert(9, 123);
    auto t = map.lookup(9);

    ASSERT_EQ(map.lookup(9), 123);
    ASSERT_EQ(map.lookup(90), ~0ul);
}

TEST(HashMap, Rehashing) {
    /*
     * the HashMap is rehashed using a copy constructor
     */
    HashMap<size_t> map(2);
    map.insert(4, 5);
    map.insert(8, 9);
    map.insert(12, 13);
    map.insert(16, 17);
    ASSERT_EQ(map.size(), 4);
    map = HashMap<size_t>(map, 4);
    ASSERT_EQ(map.lookup(4), 5);
    ASSERT_EQ(map.lookup(8), 9);
    ASSERT_EQ(map.lookup(12), 13);
    ASSERT_EQ(map.lookup(16), 17);
    ASSERT_EQ(map.size(), 4);
}

TEST(HashMap, Removal){
    HashMap<size_t> map(3);
    const size_t n = 100ul;
    for (size_t i = 0ul; i<n; ++i) {
        map.insert(i * i, i);
    }
    ASSERT_EQ(map.size(), n);
    for (size_t i = 0ul; i<n; ++i) {
        map.remove(i * i);
        ASSERT_EQ(map.size(), n-i-1);
    }
}
#endif