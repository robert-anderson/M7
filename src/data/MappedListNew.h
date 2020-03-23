//
// Created by rja on 12/03/2020.
//

#ifndef M7_MAPPEDLISTNEW_H
#define M7_MAPPEDLISTNEW_H

#include "ListNew.h"
#include "MutexVector.h"
#include "Field.h"
#include "src/hash/SafeHashMap.h"

template<typename Hashable_T>
class ListSafeHashMap;

template<typename Field_T>
class MappedListNew : public ListNew {
protected:
    ListSafeHashMap<typename Field_T::Element> m_map;
    const Field_T &m_key_field;

public:

    MappedListNew(size_t nbucket, const Field_T &key_field, size_t nsegment=1) :
            ListNew(nsegment), m_key_field(key_field), m_map(*this, nbucket) {
        assert(key_field->m_nelement==1);
    }

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
        size_t irow = ListNew::push();
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

    bool row_empty(const size_t &irow) const {
        return view<T>(irow).is_zero();
    }
};

template<typename T>
struct ListSafeHashMap : public SafeHashMap<T> {
    MappedListNew<T> &m_list;

    ListSafeHashMap(MappedListNew<T> &list, const size_t &nbucket) :
            SafeHashMap<T>(nbucket), m_list(list) {}

    T get_key(const size_t &key_index) const override {
        return *m_list.m_key_field(key_index);
    }

    void get_key(const size_t &key_index, T &key) override {
        key = *m_list.m_key_field(key_index) = key;
    }

    void set_key(const size_t &key_index, const T &key) override {
        *m_list.m_key_field(key_index) = key;
    }
};

#endif //M7_MAPPEDLISTNEW_H
