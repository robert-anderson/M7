//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include "Table.h"
#include "BitfieldHasher.h"
#include "DeterminantHasher.h"


template<typename T>
struct MapType {
    typedef std::unordered_map<T, size_t> type;
};

template<>
struct MapType<Determinant> {
    typedef std::unordered_map<Determinant, size_t, DeterminantHasher> type;
};

template<>
struct MapType<BitfieldNew> {
    typedef std::unordered_map<BitfieldNew, size_t, BitfieldHasher> type;
};


template<typename T>
class MappedTable : public Table {
    std::vector<typename MapType<T>::type> m_map;
    const size_t m_ihashentry;

public:
    MappedTable(Specification spec, size_t nrow_initial, size_t m_nsegment = 1,
                float nrow_growth_factor = 2.0, size_t nrow_mutex_blocks = 0, size_t ihashentry = 0) :
            Table(spec, nrow_initial, m_nsegment, nrow_growth_factor, nrow_mutex_blocks),
            m_map(m_nsegment), m_ihashentry(ihashentry) {}

    bool row_filled(const size_t &isegment, const size_t &irow) {
        return !view<T>(isegment, irow, m_ihashentry).is_zero();
    }

    size_t lookup(const T &key, const size_t &isegment) const {
        for (auto i=m_map[isegment].find(key); i!=m_map[isegment].end(); ++i){
            if (i->first==key) return i->second;
        }
        return ~0ul;
    }

    void insert(const T &key, size_t isegment, size_t irow) {
        m_map[isegment].emplace(key, irow);
        view<Determinant>(isegment, irow, m_ihashentry) = key;
    }

    void safe_insert(const T &key, size_t isegment, size_t irow) {
        m_segment_mutex.acquire_lock(isegment);
        insert(key, isegment, irow);
        m_segment_mutex.release_lock(isegment);
    }

    size_t push(const size_t &isegment, const T &key) {
        auto irow = lookup(key, isegment);
        if (irow == ~0ul) {
            // the key was not found in the hash map, so push a new row
            irow = Table::push(isegment, 1);
            insert(key, isegment, irow);
        }
        return irow;
    }

    size_t safe_push(const size_t &isegment, const T &key) {
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

};


#endif //M7_MAPPEDTABLE_H
