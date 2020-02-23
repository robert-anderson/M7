//
// Created by Robert John Anderson on 2020-02-19.
//

#ifndef M7_SAFEHASHMAP_H
#define M7_SAFEHASHMAP_H

#include "HashMap.h"
#include "MutexVector.h"


template<typename T>
class SafeHashMap : public HashMap<T> {
    using HashMap<T>::bucket;
    mutable MutexVector m_bucket_mutex;
public:
    SafeHashMap(const size_t &nbucket) :
        HashMap<T>(nbucket), m_bucket_mutex(nbucket) {}

    SafeHashMap(const HashMap<T> &old, const size_t &nbucket) :
        HashMap<T>(old, nbucket), m_bucket_mutex(nbucket) {}

    Mutex get_mutex(const size_t &ibucket) const {
        return m_bucket_mutex.get(ibucket);
    }

    Mutex find_mutex(const T &key) const {
        return m_bucket_mutex.get(bucket(key));
    }

    std::pair<T, size_t> *lookup_pair(Mutex mutex, const T &key) const {
        return HashMap<T>::lookup_pair(mutex.index(), key);
    }

    std::pair<T, size_t> *lookup_pair(const size_t &ibucket, const T &key) const {
        auto mutex = get_mutex(ibucket);
        return lookup_pair(mutex, key);
    }

    std::pair<T, size_t> *lookup_pair(const T &key) const {
        auto mutex = get_mutex(bucket(key));
        return lookup_pair(mutex, key);
    }

    size_t lookup(Mutex mutex, const T &key) const {
        auto tmp = lookup_pair(mutex, key);
        if (tmp) return tmp->second;
        else return ~0ul;
    }

    size_t lookup(const size_t &ibucket, const T &key) const {
        auto mutex = get_mutex(ibucket);
        return lookup(mutex, key);
    }

    size_t lookup(const T &key) const {
        auto mutex = get_mutex(bucket(key));
        return lookup(mutex, key);
    }

    void insert(Mutex mutex, const T &key, const size_t &value) {
        HashMap<T>::insert(mutex.index(), key, value);
    }

    virtual void insert(const T &key, const size_t &value) {
        auto mutex = find_mutex(key);
        return insert(mutex, key, value);
    }

    void replace(Mutex mutex, const T &key, const size_t &value) {
        auto pair = lookup_pair(mutex, key);
        if (pair) pair->second = value;
        else insert(mutex, key, value);
    }

    size_t remove(Mutex mutex, const T &key) {
        return HashMap<T>::remove(mutex.index(), key);
    }

    virtual size_t remove(const T &key) {
        auto mutex = find_mutex(key);
        return remove(mutex, key);
    }

    virtual void replace(const T &key, const size_t &value) {
        auto mutex = find_mutex(key);
        replace(mutex, key, value);
    }
};

#endif //M7_SAFEHASHMAP_H
