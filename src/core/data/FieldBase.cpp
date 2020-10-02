//
// Created by rja on 02/10/2020.
//

#include "FieldBase.h"

FieldBase::FieldBase(Table* table, size_t element_size, size_t nelement, const std::type_info& type_info) :
        m_table(table), m_element_size(element_size), m_nelement(nelement),
        m_size(m_element_size*m_nelement), m_type_info(type_info){}

bool FieldBase::is_same_type_as(const FieldBase &other) const {
    return other.m_type_info == m_type_info;
}

char *FieldBase::begin(const size_t &irow) const {
    return m_table->row_ptr(irow)+m_offset;
}
