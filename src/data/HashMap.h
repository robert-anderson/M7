//
// Created by Robert John Anderson on 2020-02-18.
//

#ifndef M7_HASHMAP_H
#define M7_HASHMAP_H

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
class HashMap {
    std::vector<std::forward_list<std::pair<T, size_t>>> m_buckets;
    typename HasherType<T>::type m_hasher;
public:
    HashMap(const size_t &nbucket) : m_buckets(nbucket){}
    HashMap(const HashMap &old, const size_t &nbucket) :
    HashMap(nbucket) {
        for (auto bucket:old.m_buckets){
            for (auto it:bucket){
                insert(it.first, it.second);
            }
        }
    }

    size_t bucket(const T &key) const {
        return m_hasher(key) % m_buckets.size();
    }

    virtual std::pair<T, size_t>* lookup_pair(const size_t &ibucket, const T &key) const {
        std::pair<T, size_t>* pair = nullptr;
        for (auto it : m_buckets[ibucket]) {
            if (it.first == key) {
                pair = &it;
                break;
            }
        }
        return pair;
    }

    virtual std::pair<T, size_t>* lookup_pair(const T &key) const {
        return lookup_pair(bucket(key), key);
    }

    virtual size_t lookup(const size_t &ibucket, const T &key) const {
        auto tmp = lookup_pair(ibucket, key);
        if (tmp) return tmp->second;
        else return ~0ul;
    }

    virtual size_t lookup(const T &key) const {
        return lookup(bucket(key), key);
    }

    virtual void insert(const size_t &ibucket, const T& key, const size_t &value) {
        /*
         * assumes the key is not already in the map
         */
        m_buckets[ibucket].emplace_front(std::pair<T, size_t>(key, value));
    }

    virtual void insert(const T& key, const size_t &value) {
        return insert(bucket(key), key, value);
    }

    void replace(const T &key, const size_t &value) {
        auto pair = lookup_pair(key);
        if (pair) pair->second = value;
        else insert(key, value);
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
            for (auto it : bucket) {
                if (it.second != other.HashMap<T>::lookup(it.first, ibucket)) return false;
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


#endif //M7_HASHMAP_H
