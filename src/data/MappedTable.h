//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include "Table.h"
#include "BitfieldHasher.h"
#include "DeterminantHasher.h"
#include "HashMap.h"

class KeyException : public std::exception {
    virtual const char *what() const throw() {
        return "Trying to access a view on a row indexed by a key which was not found in the hashmap";
    }
};

static const KeyException key_exception;

template<typename T>
class MappedTable : public Table {
protected:
    std::vector<HashMap<T>> m_map;
    const size_t m_ihashentry;
public:
    /*
    MappedTable(Specification spec, size_t nrow_initial, size_t m_nsegment = 1,
                float nrow_growth_factor = 2.0, size_t nrow_mutex_blocks = 0, size_t ihashentry = 0) :
            Table(spec, nrow_initial, m_nsegment, nrow_growth_factor, nrow_mutex_blocks),
            m_map(m_nsegment, {nrow_initial}), m_ihashentry(ihashentry) {}

    bool row_filled(const size_t &isegment, const size_t &irow) {
        return !view<T>(isegment, irow, m_ihashentry).is_zero();
    }

    size_t lookup(const T &key, const size_t &isegment) const {
        return m_map[isegment].lookup(key);
    }

    void insert(const T &key, size_t isegment, size_t irow) {
        m_map[isegment].insert(key, irow);
        view<T>(isegment, irow, m_ihashentry) = key;
    }

    virtual size_t push(const size_t &isegment, const T &key) {
        auto irow = lookup(key, isegment);
        if (irow == ~0ul) {
            // the key was not found in the hash map, so push a new row
            irow = Table::push(isegment, 1);
            insert(key, isegment, irow);
        }
        return irow;
    }
            */

    /*
    size_t safe_push(const size_t &isegment, const T &key) {
        auto mutex = m_segment_mutex.get(isegment);
        auto irow = lookup(key, isegment);
        if (irow == ~0ul) {
            // the key was not found in the hash map, so push a new row
            m_segment_mutex.acquire_lock(isegment);
            irow = Table::push(isegment, 1);
            insert(key, isegment, irow);
            m_segment_mutex.release_lock(isegment);
        }
        return irow;
    }

    template<typename data_T>
    view_t <data_T> lookup_view(const size_t &isegment, const T &key, const size_t &ientry) const {
        auto irow = lookup(key, isegment);
        if (irow == ~0ul) throw KeyException();
        return Table::view<data_T>(isegment, irow, ientry);
    }

    template<typename data_T>
    view_t <data_T> lookup_view(const size_t &isegment, const T &key) const {
        return lookup_view<data_T>(isegment, key, 0);
    }

    template<typename data_T>
    view_t <data_T> lookup_view(const T &key) const {
        return lookup_view<data_T>(0, key);
    }

    template<typename data_T>
    view_t <data_T> push_view(const size_t &isegment, const T &key, const size_t &ientry) {
        auto irow = push(isegment, key);
        return Table::view<data_T>(isegment, irow, ientry);
    }

    template<typename data_T>
    view_t <data_T> push_view(const size_t &isegment, const T &key) {
        return push_view<data_T>(isegment, key, 0);
    }

    template<typename data_T>
    view_t <data_T> push_view(const T &key) {
        return push_view<data_T>(0, key, 0);
    }

    template<typename data_T>
    view_t <data_T> safe_push_view(const size_t &isegment, const T &key, const size_t &ientry) {
        auto irow = safe_push(isegment, key);
        return Table::view<data_T>(isegment, irow, ientry);
    }

    template<typename data_T>
    view_t <data_T> safe_push_view(const size_t &isegment, const T &key) {
        return safe_push_view<data_T>(isegment, key, 0);
    }

    template<typename data_T>
    view_t <data_T> safe_push_view(const T &key) {
        return safe_push_view<data_T>(0, key, 0);
    }
*/
};


#endif //M7_MAPPEDTABLE_H
