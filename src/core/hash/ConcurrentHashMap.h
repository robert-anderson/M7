//
// Created by Robert John Anderson on 2020-08-01.
//

#ifndef M7_CONCURRENTHASHMAP_H
#define M7_CONCURRENTHASHMAP_H

#include "src/core/util/defs.h"
#include "ConcurrentLinkedList.h"
#include <vector>
#include <functional>

template<typename key_T>
struct ConcurrentHashMap;

template<typename key_T>
struct ConcurrentHashMapBucket : public ConcurrentLinkedList<size_t> {
    const ConcurrentHashMap<key_T> &m_map;
    ConcurrentLinkedList<ConcurrentLinkedListNode<size_t> *> m_tombstone_prevs;

    ConcurrentHashMapBucket(const ConcurrentHashMap<key_T> &map):m_map(map){}

    void mark_for_delete_after(ConcurrentLinkedListNode<size_t> *prev) {
        m_tombstone_prevs.append(prev);
        prev->m_next->m_payload = ~0ul;
    }

    void mark_for_delete(const size_t& index) {
        for (ConcurrentLinkedListNode<size_t> *prev = this; prev->m_next; prev = prev->m_next) {
            if (prev->m_next->m_payload==index) mark_for_delete_after(prev);
        }
    }

    void clear_tombstones() {
        for (auto prev = m_tombstone_prevs.m_next; prev; prev = prev->m_next) {
            m_map.delete_key(prev->m_payload);
            delete_after(prev->m_payload);
        }
    }

    ConcurrentLinkedListNode<size_t> *lookup_prev(const key_T &key) {
        for (ConcurrentLinkedListNode<size_t> *prev = this; prev->m_next; prev = prev->m_next) {
            if (m_map.get_key(prev->m_next->m_payload) == key) return prev;
        }
        return nullptr;
    }

    size_t lookup_index(const key_T &key) {
        auto ptr = lookup_prev(key);
        return ptr ? ptr->m_next->m_payload : ~0ul;
    }

};

template<typename key_T>
struct ConcurrentHashMap {
    size_t m_nbucket;
    std::vector<ConcurrentHashMapBucket<key_T>> m_buckets;
    std::function<void(const size_t &index)> m_delete_callback = [](const size_t &index) {};

    explicit ConcurrentHashMap(size_t nbucket) : m_nbucket(nbucket),
    m_buckets(m_nbucket, *this) {}

    void mark_for_delete(const size_t& index) {
        m_buckets[get_key(index).hash() % m_nbucket].mark_for_delete(index);
    }

    void clear_tombstones() {
        // pragma omp parallel for
        for (auto bucket : m_buckets) bucket.clear_tombstones();
    }

    ConcurrentLinkedListNode<size_t> *lookup_prev(const key_T &key) {
        return m_buckets[key.hash() % m_nbucket].lookup_prev(key);
    }

    size_t lookup_index(const key_T &key) {
        return m_buckets[key.hash() % m_nbucket].lookup_index(key);
    }

    ConcurrentLinkedListNode<size_t> *insert(const key_T &key, const size_t& index) {
        // key should not already exist in bucket
        ASSERT(!lookup_prev(key));
        return m_buckets[key.hash() % m_nbucket].append(index);
    }

    virtual key_T get_key(const size_t &index) const = 0;

    virtual void set_key(const size_t &index, const key_T &key) = 0;

};


#endif //M7_CONCURRENTHASHMAP_H
