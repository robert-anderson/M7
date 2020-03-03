//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_MAPPEDLIST_H
#define M7_MAPPEDLIST_H

#include "List.h"
#include "src/hash/SafeHashMap.h"

template<typename T>
class ListSafeHashMap;

template<typename T>
class MappedList : public List {
protected:
    const size_t m_key_entry;
    ListSafeHashMap<T> m_map;

public:
    MappedList(Specification spec, size_t nrow, size_t key_entry) :
            List(spec, nrow), m_key_entry(key_entry), m_map(*this, nrow) {}

    MappedList(Specification spec, size_t nrow, size_t key_entry, defs::data_t *data_external) :
            List(spec, nrow, data_external), m_key_entry(key_entry), m_map(*this, nrow) {}

    Mutex get_mutex(const size_t &ibucket) {
        return m_map.get_mutex(ibucket);
    }

    Mutex key_mutex(const T &key) {
        return m_map.key_mutex(key);
    }

    size_t lookup(const T &key) {
        auto mutex = key_mutex(key);
        return m_map.lookup(mutex, key);
    }

    size_t lookup(Mutex &mutex, const T &key) {
        return m_map.lookup(mutex, key);
    }

    virtual size_t push(Mutex &mutex, const T &key) {
        size_t irow = List::push();
        m_map.insert(mutex, key, irow);
        return irow;
    }

    virtual size_t push(const T &key) {
        auto mutex = key_mutex(key);
        return push(mutex, key);
    }

    size_t lookup_push(Mutex &mutex, const T &key) {
        size_t irow;
        irow = m_map.lookup(mutex, key);
        if (irow != ~0ul) return irow;
        else return push(mutex, key);
    }

    size_t lookup_push(const T &key) {
        auto mutex = key_mutex(key);
        return lookup_push(mutex, key);
    }

    const size_t &key_entry() const {
        return m_key_entry;
    }

    bool row_empty(const size_t &irow) const {
        return view<T>(irow).is_zero();
    }
};

template<typename T>
struct ListSafeHashMap : public SafeHashMap<T> {
    MappedList<T> &m_list;

    ListSafeHashMap(MappedList<T> &list, const size_t &nbucket) :
            SafeHashMap<T>(nbucket), m_list(list) {}

    T get_key(const size_t &key_index) const override {
        return m_list.template view<T>(key_index, m_list.key_entry());
    }

    void set_key(const size_t &key_index, const T &key) override {
        m_list.template view<T>(key_index, m_list.key_entry()) = key;
    }
};


#endif //M7_MAPPEDLIST_H
