//
// Created by Robert John Anderson on 2020-02-18.
//


#include <gtest/gtest.h>
#include <src/hash/SafeHashMap.h>

#if 0
TEST(SafeHashMap, WriteThreadSafety) {
    const size_t n = 10000;
    auto f = [n]() {
        SafeHashMap<size_t> map(n, KeyMap<size_t>());
#pragma omp parallel for default(none) shared(map)
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
    const size_t n = 1000;
    auto f = [n]() {
        SafeHashMap<size_t> map(n);
        SafeHashMap<size_t> read_map(n);
#pragma omp parallel for
        for (size_t i = 0ul; i < n; ++i) {
            map.insert(i, 2 * i);
        }
#pragma omp parallel for
        for (size_t i = 0ul; i < n; ++i) {
            if (i % 2) read_map.insert(i, i);
        }
        return read_map;
    };
    ASSERT_EQ(f(), f());
}

TEST(SafeHashMap, Removal) {
    SafeHashMap<size_t> map(3);
    const size_t n = 100ul;
#pragma omp parallel for
    for (size_t i = 0ul; i < n; ++i) {
        map.insert(i * i, i);
    }
    ASSERT_EQ(map.size(), n);
#pragma omp parallel for
    for (size_t i = 0ul; i < n; ++i) {
        auto mutex = map.key_mutex(i*i);
        map.remove(mutex, i);
    }
    ASSERT_EQ(map.size(), 0);
}
#endif
