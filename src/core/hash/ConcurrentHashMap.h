//
// Created by Robert John Anderson on 2020-08-01.
//

#ifndef M7_CONCURRENTHASHMAP_H
#define M7_CONCURRENTHASHMAP_H

#include "src/core/util/defs.h"
#include "ConcurrentLinkedList.h"
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <src/core/util/utils.h>

template<typename key_T>
struct ConcurrentHashMap;

template<typename key_T>
struct ConcurrentHashMapBucket : public ConcurrentLinkedList<size_t> {
#if 0
    ConcurrentLinkedList<ConcurrentLinkedListNode<size_t> *> m_tombstone_prevs;
#endif

    void mark_for_delete_after(ConcurrentLinkedListNode<size_t> *prev) {
#if 0
        m_tombstone_prevs.append(prev);
#endif
        prev->m_next->m_payload = ~0ul;
    }

    void mark_for_delete(const size_t &index) {
        for (ConcurrentLinkedListNode<size_t> *prev = this; prev->m_next; prev = prev->m_next) {
            if (prev->m_next->m_payload == index) {
#if 0
                mark_for_delete_after(prev);
#endif
                delete_after(prev);
                return;
            }
        }
        throw std::runtime_error("Index not found in hashmap bucket");
    }

#if 0
    void clear_tombstones() {
        for (auto prev = m_tombstone_prevs.m_next; prev; prev = prev->m_next) {
            delete_after(prev->m_payload);
        }
    }

    size_t ntombstone() const{
        return m_tombstone_prevs.size();
    }
#endif
};

template<typename key_T>
struct ConcurrentHashMap {
    size_t m_nbucket;
    std::vector<ConcurrentHashMapBucket<key_T>> m_buckets;

    explicit ConcurrentHashMap(size_t nbucket) : m_nbucket(nbucket), m_buckets(m_nbucket) {}

    void mark_for_delete(const size_t &index) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "MARKING HASH TABLE ENTRY FOR DELETION" << std::endl;
        std::cout << consts::verb << "index:            " << index << std::endl;
#endif
        auto key = get_key(index);
        m_buckets[hash(key) % m_nbucket].mark_for_delete(index);
        zero_key(index);
    }

#if 0
    void clear_tombstones() {
        // pragma omp parallel for
        for (auto& bucket : m_buckets) bucket.clear_tombstones();
    }
#endif

    ConcurrentLinkedListNode<size_t> *lookup_prev(const key_T &key) {
        auto &bucket = m_buckets[hash(key) % m_nbucket];
        for (ConcurrentLinkedListNode<size_t> *prev = &bucket; prev->m_next; prev = prev->m_next) {
            if (get_key(prev->m_next->m_payload) == key) return prev;
        }
        return nullptr;
    }

    size_t lookup_index(const key_T &key) {
        auto ptr = lookup_prev(key);
        return ptr ? ptr->m_next->m_payload : ~0ul;
    }

    ConcurrentLinkedListNode<size_t> *insert(const key_T &key, const size_t &index) {
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "INSERTING HASH TABLE ENTRY" << std::endl;
        std::cout << consts::verb << "index:            " << index << std::endl;
#endif
        // key should not already exist in bucket
        ASSERT(!lookup_prev(key));
        set_key(index, key);
        auto &bucket = m_buckets[hash(key) % m_nbucket];
        auto prev = bucket.append(index);
        ASSERT(!bucket.is_empty())
        return prev;
    }

    std::vector<size_t> to_vector() const {
        std::vector<size_t> out;
        for (const auto &bucket : m_buckets) {
            if (!bucket.is_empty()) {
                auto tmp = bucket.to_vector();
                out.insert(out.end(), tmp.begin(), tmp.end());
            }
        }
        return out;
    }

    size_t size() const {
        return std::accumulate(
            m_buckets.begin(), m_buckets.end(), 0,
            [](size_t acc, const ConcurrentHashMapBucket<key_T> &bucket) {
                return acc + bucket.size();
            }
        );
    }

#if 0
    size_t ntombstone() const {
        return std::accumulate(
                m_buckets.begin(), m_buckets.end(), 0,
                [](size_t acc, const ConcurrentHashMapBucket<key_T> &bucket) {
                    return acc + bucket.ntombstone();
                }
        );
    }
#endif

    bool is_empty() const {
        return std::all_of(m_buckets.begin(), m_buckets.end(),
                           [](const ConcurrentHashMapBucket<key_T> &bucket) { return bucket.is_empty(); });
    }

    virtual size_t hash(const key_T &key) const = 0;

    virtual key_T get_key(const size_t &index) const = 0;

    virtual void set_key(const size_t &index, const key_T &key) = 0;

    virtual void zero_key(const size_t &index) = 0;
};


#endif //M7_CONCURRENTHASHMAP_H
