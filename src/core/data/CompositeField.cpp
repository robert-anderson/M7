//
// Created by RJA on 01/11/2020.
//

#include "CompositeField.h"

CompositeField::CompositeField(TableX *table) : m_table(table) {}

size_t CompositeField::add_field(const NdFieldBaseX *field) {
    FieldGroup::add_field(field);
    auto offset = m_table->add_field(field);
    return offset;
}

const char *CompositeField::raw_ptr(const size_t &icomponent, const size_t &irow, const size_t &ielement) const {
    ASSERT(icomponent<nfield())
    return m_fields[icomponent]->raw_ptr(irow, ielement);
}

const size_t &CompositeField::size(const size_t &icomponent) const {
    return m_fields[icomponent]->m_size;
}

const size_t &CompositeField::offset(const size_t &icomponent) const {
    return m_fields[icomponent]->m_offset;
}
