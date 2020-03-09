//
// Created by rja on 05/03/2020.
//

#ifndef M7_FIELDSET_H
#define M7_FIELDSET_H

#include <src/defs.h>
#include <iostream>
#include <assert.h>
#include <src/multidim/Indexer.h>

struct FieldBase;

template <typename T>
struct Field;

/*
 * FieldSets both define the layout of data stored in a defs::data_t buffer,
 * and provides a namespace in which the stored data can be accessed by
 * variables with the Field type.
 *
 * Specification of the data layout is done by building up the member fields
 * in the list initialization of the FieldSet struct. Consecutively added 
 * fields of the same data type are allowed to share a defs::data_t word,
 * with the exception of boolean fields of length>1 (aka bitfields). Bitfields
 * are always allotted a whole number of data words.
 *
 * e.g. A: 3 floats
 *      B: 5 chars
 *      C: 6 chars
 *      D: 28 bools
 *      E: 1 bool
 *      F: 1 bool
 *      G: 1 bool
 *
 * |A0A0A0A0|A0A0A0A0|A0A0A0A0|A0A0A0A0|A1A1A1A1|A1A1A1A1|A1A1A1A1|A1A1A1A1|
 * |A2A2A2A2|A2A2A2A2|A2A2A2A2|A2A2A2A2|________|________|________|________|
 * |C0C0C0C0|C1C1C1C1|C2C2C2C2|C3C3C3C3|C4C4C4C4|C5C5C5C5|________|________|
 * |DDDDDDDD|DDDDDDDD|DDDDDDDD|DDDD____|________|________|________|________|
 * |EFG_____|________|________|________|________|________|________|________|
 *
 */

struct FieldSet {
    defs::data_t *m_buffer = nullptr;
    size_t m_length = 0ul;

    template<typename T>
    size_t add_field(Field<T> *field) {
    }

    /*
    template<typename T>
    size_t add_field(Field<T> *field) {
        auto advance = [&]() {
            auto tmp = m_nbit % size_in_bits<defs::data_t>();
            if (tmp) m_nbit += size_in_bits<defs::data_t>() - tmp;
        };

        size_t this_type = typeid(T).hash_code();
        if (m_nbit > 0) {
            if (this_type != m_last_type) advance();
        }

        auto offset = m_nbit / size_in_bits<T>();
        m_nbit += field->m_length * size_in_bits<T>();
        m_length = m_nbit / size_in_bits<defs::data_t>();
        if (m_nbit % size_in_bits<defs::data_t>()) m_length++;
        m_last_type = this_type;
        m_fields.push_back(field);
        return offset;
    }
     */

    std::vector<FieldBase *> m_fields;

    FieldSet(defs::data_t *buffer) : m_buffer(buffer) {}

    defs::inds field_lengths();

    defs::inds field_offsets();

    void print(size_t irow);
};


#endif //M7_FIELDSET_H
