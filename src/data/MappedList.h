//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_MAPPEDLIST_H
#define M7_MAPPEDLIST_H

#include "List.h"
#include "SafeHashMap.h"

template<typename T>
class MappedList : public List {
protected:
    SafeHashMap<T> m_map;
    const size_t m_key_entry;
public:
    MappedList(Specification spec, size_t nrow, size_t key_entry) :
    List(spec, nrow), m_map(nrow), m_key_entry(key_entry) {}
    MappedList(Specification spec, size_t nrow, size_t key_entry, defs::data_t *data_external) :
    List(spec, nrow, data_external), m_map(nrow), m_key_entry(key_entry) {}

    Mutex get_mutex(const size_t &ibucket){
        return m_map.get_mutex(ibucket);
    }

    Mutex find_mutex(const T& key){
        return m_map.find_mutex(key);
    }

    size_t lookup(const T &key){
        auto mutex = find_mutex(key);
        return m_map.lookup(mutex, key);
    }

    size_t lookup(Mutex &mutex, const T &key){
        return m_map.lookup(mutex, key);
    }

    size_t push(Mutex &mutex, const T &key){
        size_t irow;
        irow = m_map.lookup(mutex, key);
        if (irow!=~0ul) return irow;
        irow = List::push();
        m_map.insert(mutex, key, irow);
        view<T>(irow, m_key_entry) = key;
        return irow;
    }

    size_t push(const T& key){
        auto mutex = find_mutex(key);
        return push(mutex, key);
    }
};


#endif //M7_MAPPEDLIST_H
