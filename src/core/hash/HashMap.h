//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef M7_HASHMAP_H
#define M7_HASHMAP_H


#include <vector>
#include <forward_list>
#include <src/core/table/Element.h>
#include <iostream>
#include "Hashing.h"

/*
 * HashMap maps a set of keys onto a set of integral indices via a hash function.
 * Collisions are resolved by iterating through a singly linked list of
 * key index, value index pairs
 *
 * "lookup" refers to searching through the entries of a bucket until a matching
 * key is found. Returns null pair or ~0ul if no match is found
 *
 * "insert" places the new key index, value index pair in the map, and places
 * the new key in the KeyMap. Assumes the key is not already in the map. It is
 * the caller's responsibility to perform a prior lookup, or in some other way
 * make sure that the key index is not already in the map.

 * "remove" finds and deletes the pair node matching a key index, and returns
 * the removed value.
 */

template<typename T>
class HashMap {
    static_assert(std::is_base_of<Element, T>::value,
        "Only types derived from Element can serve as the key type for a HashMap");
protected:
    std::vector<std::forward_list<size_t>> m_buckets;
public:
    HashMap(const size_t &nbucket) : m_buckets(nbucket) {}

    HashMap(const HashMap &old, const size_t &nbucket) : HashMap(nbucket) {
        size_t ibucket = 0ul;
        for (auto old_bucket:old.m_buckets) {
            for (auto it:old_bucket) {
                m_buckets[bucket(old.get_key(it))].emplace_front(it);
            }
            ++ibucket;
        }
    }

    virtual T get_key(const size_t &key_index) const = 0;

    virtual void set_key(const size_t &key_index, const T &key) = 0;

    size_t bucket(const T& key) const {
        return key.hash() % m_buckets.size();
    }

    virtual size_t lookup(const size_t &ibucket, const T &key) const {
        auto &list = m_buckets[ibucket];
        for (auto &it : list) {
            if (get_key(it) == key) return it;
        }
        return ~0ul;
    }

    virtual size_t lookup(const T &key) const {
        return lookup(bucket(key), key);
    }

    virtual void insert(const size_t &ibucket, const T &key, const size_t &key_index) {
        m_buckets[ibucket].emplace_front(key_index);
        set_key(key_index, key);
    }

    virtual void insert(const T &key, const size_t &key_index) {
        insert(bucket(key), key, key_index);
    }

    virtual size_t remove(const size_t &ibucket, const size_t &key_index) {
        auto &list = m_buckets[ibucket];
        auto prev = list.begin();
        for (auto it = list.begin(); it != list.end(); it++) {
            auto t = *it;
            if (*it == key_index) {
                if (prev == it) list.pop_front();
                else list.erase_after(prev);
                return key_index;
            }
            prev = it;
        }
        return ~0ul;
    }

    virtual size_t remove(const T &key) {
        auto &list = m_buckets[bucket(key)];
        auto prev = list.begin();
        for (auto it = list.begin(); it != list.end(); it++) {
            if (get_key(*it) == key) {
                if (prev == it) list.pop_front();
                else list.erase_after(prev);
                return *it;
            }
            prev = it;
        }
        return ~0ul;
    }

    size_t size() const {
        size_t tmp = 0;
        for (auto bucket : m_buckets) {
            tmp += std::distance(bucket.begin(), bucket.end());
        }
        return tmp;
    }

    bool operator==(const HashMap &other) const {
        if (size() != other.size()) return false;
        auto ibucket = 0ul;
        for (auto bucket : m_buckets) {
            for (auto &it : bucket) {
                if (it != other.lookup(ibucket, get_key(it))) return false;
            }
            ibucket++;
        }
        return true;
    }

    void print() const {
        size_t ibucket=0ul;
        for (auto &bucket : m_buckets) {
            std::cout << std::endl << "Bucket " << ibucket << std::endl;
            for (auto &it : bucket) {
                auto tmp = get_key(it);
                auto str = get_key(it).to_string();
                std::cout << get_key(it).to_string() << "->" << it << std::endl;
            }
            ibucket++;
        }
    }
};


#endif //M7_HASHMAP_H
