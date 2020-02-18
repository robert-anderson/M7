//
// Created by Robert John Anderson on 2020-02-18.
//

#ifndef M7_SAFEHASHMAP_H
#define M7_SAFEHASHMAP_H

#include <iostream>
#include <functional>
#include <forward_list>
#include <src/fermion/Determinant.h>
#include "DeterminantHasher.h"
#include "BitfieldHasher.h"
#include "MutexVector.h"

template<typename T>
struct HasherType {
    typedef std::hash<T> type;
};

template<>
struct HasherType<Determinant> {
    typedef DeterminantHasher type;
};

template<>
struct HasherType<BitfieldNew> {
    typedef BitfieldHasher type;
};

template<typename T>
class SafeHashMap {
    std::vector<std::forward_list<std::pair<T, size_t>>> m_buckets;
    typename HasherType<T>::type m_hasher;
    MutexVector m_bucket_mutex;
public:
    SafeHashMap(const size_t &nbucket) : m_buckets(nbucket), m_bucket_mutex(nbucket) {}
    SafeHashMap(const SafeHashMap &old, const size_t &nbucket) :
    SafeHashMap(nbucket) {
        for (auto bucket:old.m_buckets){
            for (auto it:bucket){
                insert(it.first, it.second);
            }
        }
    }

    size_t bucket(const T &key) const {
        return m_hasher(key) % m_buckets.size();
    }

private:
    std::pair<T, size_t>* no_mutex_lookup_pair(const T &key, const size_t &ibucket) const {
        std::pair<T, size_t>* pair = nullptr;
        for (auto it : m_buckets[ibucket]) {
            if (it.first == key) {
                pair = &it;
                break;
            }
        }
        return pair;
    }

    size_t no_mutex_lookup(const T &key, const size_t &ibucket) const {
        auto tmp = no_mutex_lookup_pair(key, ibucket);
        if (tmp) return tmp->second;
        else return ~0ul;
    }

    void no_mutex_insert(const T &key, const size_t &value, const size_t &ibucket) {
        /*
         * assumes the key is not already in the map
         */
        m_buckets[ibucket].emplace_front(std::pair<T, size_t>(key, value));
    }

public:

    std::pair<T, size_t>* lookup_pair(const T &key) {
        std::pair<T, size_t>* pair = nullptr;
        size_t ibucket = bucket(key);
        m_bucket_mutex.acquire_lock(ibucket);
        for (auto it : m_buckets[ibucket]) {
            if (it.first == key) {
                pair = &it;
                break;
            }
        }
        m_bucket_mutex.release_lock(ibucket);
        return pair;
    }

    size_t lookup(const T &key) {
        auto tmp = lookup_pair(key);
        if (tmp) return tmp->second;
        else return ~0ul;
    }

    void insert(const T &key, const size_t &value) {
        size_t ibucket = bucket(key);
        m_bucket_mutex.acquire_lock(ibucket);
        no_mutex_insert(key, value, ibucket);
        m_bucket_mutex.release_lock(ibucket);
    }

    void replace(const T &key, const size_t &value) {
        size_t ibucket = bucket(key);
        m_bucket_mutex.acquire_lock(ibucket);
        auto pair = no_mutex_lookup_pair(key, ibucket);
        if (pair) pair->second = value;
        else no_mutex_insert(key, value, ibucket);
        m_bucket_mutex.release_lock(ibucket);
    }


    size_t size() const {
        /*
         * not thread safe
         */
        size_t tmp = 0;
        for (auto bucket : m_buckets) {
            tmp += std::distance(bucket.begin(), bucket.end());
        }
        return tmp;
    }

    bool operator==(const SafeHashMap<T> &other) const {
        if (size() != other.size()) return false;
        auto ibucket = 0ul;
        for (auto bucket : m_buckets) {
            for (auto it : bucket) {
                if (it.second != other.no_mutex_lookup(it.first, ibucket)) return false;
            }
            ibucket++;
        }
        return true;
    }

    void print() const {
        for (auto bucket : m_buckets) {
            for (auto it : bucket) {
                std::cout << it.first << "->" << it.second << std::endl;
            }
        }
    }
};


#endif //M7_SAFEHASHMAP_H
