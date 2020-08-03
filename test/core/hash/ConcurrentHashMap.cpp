//
// Created by Robert John Anderson on 2020-07-30.
//

#include "gtest/gtest.h"
#include <omp.h>
#include <src/core/hash/ConcurrentHashMap.h>
#include <src/core/util/utils.h>


TEST(ConcurrentHashMap, Test) {
    struct TestHashMap : public ConcurrentHashMap<size_t> {
        std::vector<size_t>& m_keys;
        TestHashMap(size_t nbucket, std::vector<size_t>& keys):
        ConcurrentHashMap<size_t>(nbucket), m_keys(keys){}
        size_t get_key(const size_t &index) const override {
            return m_keys[index];
        }

        void set_key(const size_t &index, const size_t &key) override {
            m_keys[index] = key;
        }

        size_t hash(const size_t &key) const override {
            return ((key*123ul)<<1ul)^124576781232133ul ;
        }
    };

    std::vector<size_t> keys(100);
    TestHashMap map(10, keys);
    ASSERT_TRUE(map.is_empty());
    ASSERT_EQ(map.size(), 0);

    map.insert(3, 0);
    ASSERT_EQ(map.size(), 1);
    map.insert(33, 1);
    ASSERT_EQ(map.size(), 2);
    map.insert(333, 6);
    ASSERT_EQ(map.size(), 3);
    map.insert(3333, 62);
    ASSERT_EQ(map.size(), 4);

    ASSERT_EQ(map.lookup_index(3), 0);
    ASSERT_EQ(map.lookup_index(33), 1);
    ASSERT_EQ(map.lookup_index(333), 6);
    ASSERT_EQ(map.lookup_index(3333), 62);

    map.mark_for_delete(6);
    map.mark_for_delete(1);
    ASSERT_EQ(map.size(), 4);
    map.clear_tombstones();
    ASSERT_EQ(map.size(), 2);

    map.insert(35, 1);
    ASSERT_EQ(map.size(), 3);
}



