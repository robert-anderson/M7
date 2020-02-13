//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include "Table.h"

template<typename T>
class MappedTable : public Table {
    //std::unordered_map<Determinant, size_t> m_map;
    std::unordered_map<T, size_t> m_map;
    const size_t m_ihashentry;

public:
    MappedTable(Specification spec, size_t ihashentry, size_t nrow_initial, size_t m_nsegment = 1,
        float nrow_growth_factor = 2.0, size_t nrow_mutex_blocks=0);

    bool row_filled(const size_t &isegment, const size_t &irow){
        return !view<T>(isegment, irow, m_ihashentry).is_zero();
    }

    size_t push(const size_t &isegment);

    size_t safe_push(const size_t &isegment);
};

template<>
class MappedTable<size_t>;

template<>
class MappedTable<BitfieldNew>;

template<>
class MappedTable<Determinant>;

#endif //M7_MAPPEDTABLE_H
