//
// Created by Robert John Anderson on 2020-02-18.
//


#include <gtest/gtest.h>
#include <src/data/SafeHashMap.h>


TEST(SafeHashMap, SimpleLookup) {
    const size_t n = 100;
    SafeHashMap<size_t> map(n);
    map.insert(9, 123);
    ASSERT_EQ(map.lookup(9), 123);
    ASSERT_EQ(map.lookup(90), ~0ul);
}

TEST(SafeHashMap, WriteThreadSafety) {
    const size_t n = 10000;
    auto f = [n]() {
        SafeHashMap<size_t> map(n);
#pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            map.insert(i, 2 * i);
        }
        return map;
    };
    ASSERT_EQ(f(), f());
}

TEST(SafeHashMap, ReadWriteThreadSafety) {
    /*
     * check the integrity of concurrent read/write on a SafeHashMap
     */
    const size_t n = 10000;
    auto f = [n]() {
        SafeHashMap<size_t> map(n);
        SafeHashMap<size_t> read_map(n);
#pragma omp parallel for
        for (size_t i = 0ul; i < n; ++i) {
            // filling half the map with integers
            map.insert(i, 2 * i);
        }
#pragma omp parallel for
        for (size_t i = 0ul; i < n; ++i) {
            if (i % 2) map.replace(2 * i, 3 * i);
            else {
                read_map.replace(map.lookup(i), i);
            }
        }
        return read_map;
    };
    ASSERT_EQ(f(), f());
}

TEST(SafeHashMap, Rehashing) {
    /*
     * the SafeHashMap is rehashed using a copy constructor
     */
    SafeHashMap<size_t> map(2);
    map.insert(4, 5);
    map.insert(8, 9);
    map.insert(12, 13);
    map.insert(16, 17);
    map = SafeHashMap<size_t>(map, 4);
    ASSERT_EQ(map.lookup(4), 5);
    ASSERT_EQ(map.lookup(8), 9);
    ASSERT_EQ(map.lookup(12), 13);
    ASSERT_EQ(map.lookup(16), 17);
}