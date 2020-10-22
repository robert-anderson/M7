//
// Created by rja on 21/10/2020.
//

#include "Field.h"
#include "Table.h"

FieldX::FieldX(TableX *table, size_t nelement, size_t element_size, std::string description) :
        m_table(table), m_nelement(nelement), m_element_size(element_size), m_size(nelement * element_size),
        m_offset(table->add_field(*this)), m_description(std::move(description)) {}

char *FieldX::begin(const size_t &irow) const {
    return m_table->begin(irow);
}

char *FieldX::raw_ptr(const size_t &irow, const size_t &iflat) const {
    return begin(irow) + iflat * m_element_size;
}