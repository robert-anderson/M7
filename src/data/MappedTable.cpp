//
// Created by Robert John Anderson on 2020-02-04.
//

#include "MappedTable.h"

template<typename T>
MappedTable<T>::MappedTable(Specification spec, size_t ihashentry, size_t nrow_initial, size_t m_nsegment,
                         float nrow_growth_factor, size_t nrow_mutex_blocks) :
        Table(spec, nrow_initial, m_nsegment, nrow_growth_factor, nrow_mutex_blocks),
        m_ihashentry(ihashentry){
    //m_map{};
}

template<typename T>
size_t MappedTable<T>::push(const size_t &isegment) {
    return Table::push(isegment, 1);
}

template<typename T>
size_t MappedTable<T>::safe_push(const size_t &isegment) {
    return Table::safe_push(isegment, 1);
}
