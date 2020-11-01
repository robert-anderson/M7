//
// Created by RJA on 01/11/2020.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include "Table.h"

struct CompositeField: FieldGroup {

    TableX* m_table;
    size_t m_element_size;
    size_t m_offset = ~0ul;
    CompositeField(TableX *table): m_table(table) {}

    size_t add_field(const NdFieldBaseX *field) override {
        FieldGroup::add_field(field);
        auto offset = m_table->add_field(field);
        if (m_fields.size()==1) {
            m_element_size = field->m_element_size;
            m_offset = offset;
        }
        else m_element_size = (offset-m_fields[0]->m_offset)+field->m_element_size;
        return offset;
    }

    const char *raw_ptr(const size_t &irow, const size_t &ielement) const {
        return m_table->begin(irow) + m_offset + ielement * m_element_size;
    }

    std::pair<const char *, size_t> raw_view(const size_t &irow, const size_t &ielement) const {
        return {raw_ptr(irow, ielement), m_element_size};
    }
};


#endif //M7_COMPOSITEFIELD_H
