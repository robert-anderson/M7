//
// Created by rja on 01/11/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include "Table.h"
#include "NdField.h"

//template<typename field_t, typename hash_fn=typename field_t::default_hash_fn>
//class MappedTable : TableX {
//
//    typedef typename std::conditional<
//            std::is_base_of<FieldBaseX, field_t>::value,
//            NdFieldX<field_t, 0ul>, field_t
//    >::type key_field_t;
//
//    const key_field_t &m_key_field;
//    std::unordered_map<typename field_t::view_t, size_t, hash_fn, NdFieldBaseX::equals_fn> m_map;
//
//    MappedTable(const key_field_t &key_field) : TableX(), m_key_field(key_field) {}
//
//};


#endif //M7_MAPPEDTABLE_H
