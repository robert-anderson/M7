//
// Created by Robert John Anderson on 2020-02-19.
//

#ifndef M7_SAFEHASHMAP_H
#define M7_SAFEHASHMAP_H

#include "HashMap.h"
#include "src/data/MutexVector.h"


template<typename T>
class SafeHashMap : public HashMap<T> {
    using HashMap<T>::m_buckets;
    using HashMap<T>::bucket;
    mutable MutexVector m_bucket_mutex;
public:
    SafeHashMap(const size_t &nbucket):
            HashMap<T>(nbucket), m_bucket_mutex(nbucket) {}

    SafeHashMap(const HashMap<T> &old, const size_t &nbucket) :
            HashMap<T>(old, nbucket), m_bucket_mutex(nbucket) {}

    Mutex get_mutex(const size_t &ibucket) const {
        return m_bucket_mutex.get(ibucket);
    }

    Mutex key_mutex(const T &key) const {
        return m_bucket_mutex.get(bucket(key));
    }

    size_t lookup(Mutex &mutex, const T &key) const {
        return HashMap<T>::lookup(mutex.index(), key);
    }

    size_t lookup(const T &key) const override {
        auto mutex = key_mutex(key);
        return lookup(mutex, key);
    }

    void insert(Mutex &mutex, const T &key, const size_t &key_index) {
        HashMap<T>::insert(mutex.index(), key, key_index);
    }

    void insert(const T &key, const size_t &key_index) override {
        auto mutex = key_mutex(key);
        insert(mutex, key, key_index);
    }

    size_t remove(Mutex &mutex, const size_t &key_index) {
        return HashMap<T>::remove(mutex.index(), key_index);
    }

};

#endif //M7_SAFEHASHMAP_H
