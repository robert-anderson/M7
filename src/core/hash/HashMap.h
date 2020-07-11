/**
 * @file
 * @author Robert John Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 * A HashMap forms the basis for associative random access of Table-derived objects.
 * A "key" is an object derived from the Element class, which maps to an index "irow"
 * via this class.
 *
 * The index irow refers to the position of the unique row containing the key identifier,
 * and is determined by computing the bucket to which the key belongs, then iterating
 * over a forward_list containing irows corresponding to all keys in the same bucket.
 *
 *
 *
 * consists of a number of "buckets" in which keys are
 */

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

    virtual T get_key(const size_t &irow) const = 0;

    virtual void set_key(const size_t &irow, const T &key) = 0;

    size_t bucket(const T &key) const {
        return key.hash() % m_buckets.size();
    }

    size_t lookup(const T &key) const {
        auto &list = m_buckets[bucket(key)];
        for (const auto &irow : list) {
            if (get_key(irow) == key) return irow;
        }
        return ~0ul;
    }

    void insert(const T &key, const size_t &irow) {
        m_buckets[bucket(key)].emplace_front(irow);
        set_key(irow, key);
    }

    size_t remove(const T &key, const size_t &irow) {
        auto &list = m_buckets[bucket(key)];
        auto prev = list.begin();
        for (auto it = list.begin(); it != list.end(); it++) {
            if (*it == irow) {
                if (prev == it) list.pop_front();
                else list.erase_after(prev);
                return irow;
            }
            prev = it;
        }
        return ~0ul;
    }

    size_t remove(const T &key) {
        auto irow = lookup(key);
        if (irow==~0ul) return irow;
        return remove(key, irow);
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
        size_t ibucket = 0ul;
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
