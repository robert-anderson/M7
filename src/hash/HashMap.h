//
// Created by Robert John Anderson on 2020-02-18.
//

#ifndef M7_HASHMAP_H
#define M7_HASHMAP_H

#include <iostream>
#include <memory>
#include <forward_list>
#include "DeterminantHasher.h"
#include "BitfieldHasher.h"

#if 0
#include "src/data/MutexVector.h"

template<typename T>
struct HasherType {
    typedef std::hash<T> type;
};

/*
template<>
struct HasherType<typename DeterminantField<>::Element> {
    typedef DeterminantHasher type;
};

template<>
struct HasherType<BitfieldNew> {
    typedef BitfieldHasher type;
};
*/

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
protected:
    std::vector<std::forward_list<size_t>> m_buckets;
    typename HasherType<T>::type m_hasher;
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
    virtual void set_key(const size_t &key_index, const T& key) = 0;


    size_t bucket(const T &key) const {
        return m_hasher(key) % m_buckets.size();
    }

    virtual size_t lookup(const size_t &ibucket, const T &key) const {
        auto &list = m_buckets[ibucket];
        for (auto &it : list) {
            auto y = get_key(it);
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

    bool operator==(const HashMap<T> &other) const {
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
        for (auto bucket : m_buckets) {
            for (auto &it : bucket) {
                std::cout << get_key(it) << "->" << it << std::endl;
            }
        }
    }

};


#endif //M7_HASHMAP_H
#endif //M7_HASHMAP_H
