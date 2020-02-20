//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_MAPPEDLIST_H
#define M7_MAPPEDLIST_H

#include "List.h"
#include "SafeHashMap.h"

template<typename T>
class MappedList : public List {
    SafeHashMap<T> m_map;
    const size_t m_key_entry;
public:
    MappedList(Specification spec, size_t nrow, size_t key_entry) : List(spec, nrow) {}
    MappedList(Specification spec, size_t nrow, size_t key_entry, defs::data_t *data_external) :
    List(spec, nrow, data_external), m_key_entry(key_entry) {}

    size_t push(const T& key){
        auto mutex = m_map.find_mutex(key);
        size_t irow;
        irow = m_map.lookup(mutex, key);
        if (irow!=~0ul) return irow;
        irow = List::push();
        m_map.insert(mutex, key, irow);
        return irow;
    }
};


#endif //M7_MAPPEDLIST_H
