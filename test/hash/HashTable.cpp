//
// Created by Robert John Anderson on 2020-02-20.
//


#include <gtest/gtest.h>
#include "src/hash/HashTable.h"

TEST(HashTable, SimpleLookup) {
    const size_t n = 100;
    HashTable<std::string, size_t> dict(n);
    ASSERT_EQ(dict.size(), 0);
    dict.set("nine", 9);
    ASSERT_EQ(dict.size(), 1);
    dict.set("ninety", 90);
    ASSERT_EQ(dict.size(), 2);
    ASSERT_EQ(dict.get("nine"), 9);
    ASSERT_EQ(dict.get("ninety"), 90);
    ASSERT_THROW(dict.get("seven"), KeyError);
    dict.remove("nine");
    ASSERT_THROW(dict.get("nine"), KeyError);
}

TEST(HashTable, Rehashing) {
    /*
     * the HashMap is rehashed using a copy constructor
     */
    HashTable<size_t, size_t> map(4);
    map.set(4, 5);
    map.set(8, 9);
    map.set(12, 13);
    map.set(16, 17);
    ASSERT_EQ(map.size(), 4);
    map = HashTable<size_t, size_t>(map, 6);
    ASSERT_EQ(map.get(4), 5);
    ASSERT_EQ(map.get(8), 9);
    ASSERT_EQ(map.get(12), 13);
    ASSERT_EQ(map.get(16), 17);
    ASSERT_EQ(map.size(), 4);
}

TEST(HashTable, Removal){
    HashTable<size_t, size_t> map(3);
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