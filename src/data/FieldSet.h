//
// Created by rja on 05/03/2020.
//

#ifndef M7_FIELDSET_H
#define M7_FIELDSET_H

#include <src/defs.h>
#include <iostream>
#include <assert.h>
#include <src/multidim/Indexer.h>
#include <typeindex>

/*
 * FieldSets both define the layout of data stored in a defs::data_t buffer,
 * and provide a namespace in which the stored data can be accessed by
 * symbols with the Field type.
 *
 * Specification of the data layout is done by building up the member fields
 * in the list initialization of the FieldSet struct. Consecutively added 
 * fields of the same data type are allowed to share a defs::data_t word,
 * with the exception of boolean fields of length >1 (aka bitfields). Bitfields
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
 *
 */

struct FieldSet {
    defs::data_t *m_buffer = nullptr;
    size_t m_length = 0ul;

    template<typename T>
    static size_t nbit_per_element() {
        return std::is_same<T, bool>::value ? 1 : sizeof(T) * 8;
    }

    struct FieldBase {
        FieldSet *m_field_set;
        // length of the indexer
        const size_t m_nelement;
        const size_t m_nbit_per_element;
        const std::type_index m_type_index;
        // first value is dataword offset, second is offset in elements from start of dataword
        const std::pair<size_t, size_t> m_offset;

        virtual std::string to_string(size_t irow) = 0;

        FieldBase(FieldSet *field_set, size_t nbit, size_t nelement, const std::type_info &type_info) :
                m_field_set(field_set), m_nbit_per_element(nbit), m_nelement(nelement),
                m_type_index(type_info), m_offset(field_set->add_field(this)) {}

    };

    std::vector<FieldBase *> m_fields;

    std::pair<size_t, size_t> add_field(FieldBase *field) {
        size_t dataword_offset = 0ul;
        size_t element_offset = 0ul;
        if (!m_fields.empty()) {
            if (field->m_type_index != m_fields.back()->m_type_index) {
                // different type to last field
                dataword_offset = m_length;
                element_offset = 0;
            } else {
                // same type as last field
                element_offset = m_fields.back()->m_offset.second + m_fields.back()->m_nelement;
                element_offset *= field->m_nbit_per_element;
                dataword_offset = m_fields.back()->m_offset.first;
                dataword_offset += element_offset / nbit_per_element<defs::data_t>();
                element_offset %= nbit_per_element<defs::data_t>();
            }
        }
        m_length = dataword_offset * nbit_per_element<defs::data_t>();
        m_length += (element_offset + field->m_nelement) * field->m_nbit_per_element;
        m_length = integer_utils::divceil(m_length, nbit_per_element<defs::data_t>());
        m_fields.push_back(field);
        return {dataword_offset, element_offset};
    }

    FieldSet(defs::data_t *buffer) : m_buffer(buffer) {}

    defs::inds field_nelements();

    std::vector<std::pair<size_t, size_t>> field_offsets();

    void print(size_t irow);
};


#endif //M7_FIELDSET_H
