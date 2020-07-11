//
// Created by Robert John Anderson on 2020-04-01.
//

#ifndef M7_MAPPEDLIST_H
#define M7_MAPPEDLIST_H

#include <src/core/hash/HashMap.h>
#include <src/core/table/NumericField.h>
#include "src/core/list/List.h"

template<typename T>
class MappedList;

template<typename T>
struct ListHashMap : public HashMap<T> {
    MappedList<T> &m_list;

    ListHashMap(MappedList<T> &list, const size_t &nbucket) :
            HashMap<T>(nbucket), m_list(list) {}

    T get_key(const size_t &key_index) const override {
        return m_list.key_field()(key_index);
    }

    void set_key(const size_t &key_index, const T &key) override {
        m_list.key_field()(key_index) = key;
    }
};

template<typename T>
class MappedList : public List {
    typedef typename T::Field_T Field_T;
protected:
    Field_T &m_key_field;
    ListHashMap<T> m_map;

public:
    MappedList(std::string name, Field_T &key_field, size_t nbucket) : List(name),
                                                     m_key_field(key_field), m_map(*this, nbucket) {}

    Field_T &key_field() const {
        return m_key_field;
    };

    size_t lookup(const T &key) {
        return m_map.lookup(key);
    }

    virtual size_t push(const T &key) {
        size_t irow = List::push();
        m_map.insert(key, irow);
        return irow;
    }

    size_t lookup_push(const T &key) {
        size_t irow;
        irow = m_map.lookup(key);
        if (irow != ~0ul) return irow;
        else return push(key);
    }


    bool row_empty(const size_t &irow) const {
        return m_key_field(irow).is_zero();
    }

    void print_map() const {
        m_map.print();
    }

    size_t expand_push(const T &key, const size_t &isegment, const size_t &nrow, double factor = 1.5) {
        // convenient but inefficient
        ASSERT(factor >= 1.0);
        if (high_water_mark(isegment) + nrow > nrow_per_segment()) {
            resize(factor * double(high_water_mark(isegment) + nrow));
        }
        return push(key);
    }

    size_t expand_push(const T &key) {
        return expand_push(key, 0, 1);
    }
};

#endif //M7_MAPPEDLIST_H
