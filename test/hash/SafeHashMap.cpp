//
// Created by Robert John Anderson on 2020-02-18.
//


#include <gtest/gtest.h>
#include <src/hash/SafeHashMap.h>

template <typename T>
class SafeHashTable : public SafeHashMap<T> {
    std::vector<T> m_keys;
public:
     SafeHashTable(const size_t &nbucket):
     SafeHashMap<T>(nbucket), m_keys(nbucket){}

    T get_key(const size_t &key_index) const override {
        return m_keys[key_index];
    }

    void set_key(const size_t &key_index, const T &key) override {
         m_keys[key_index] = key;
    }
};

TEST(SafeHashMap, WriteThreadSafety) {
    const size_t n = 10000;
    auto f = [n]() {
        SafeHashTable<size_t> map(n);
#pragma omp parallel for default(none) shared(map)
        for (size_t i = 0; i < n; ++i) {
            map.insert(i*3, i);
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
        SafeHashTable<size_t> map(n);
        SafeHashTable<size_t> read_map(n);
#pragma omp parallel for default(none) shared(map)
        for (size_t i = 0ul; i < n; ++i) {
            map.insert(2*i, i);
        }
#pragma omp parallel for default(none) shared(read_map)
        for (size_t i = 0ul; i < n; ++i) {
            if (i % 2) read_map.insert(i, i);
        }
        return read_map;
    };
    ASSERT_EQ(f(), f());
}


TEST(SafeHashMap, Removal) {
    const size_t n = 100ul;
    SafeHashTable<size_t> map(n);
#pragma omp parallel for default(none) shared(map)
    for (size_t i = 0ul; i < n; ++i) {
        map.insert(i * i, i);
    }
    ASSERT_EQ(map.size(), n);
#pragma omp parallel for default(none) shared(map)
    for (size_t i = 0ul; i < n; ++i) {
        auto mutex = map.key_mutex(i*i);
        map.remove(mutex, i);
    }
    ASSERT_EQ(map.size(), 0);
}