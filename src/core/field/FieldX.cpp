//
// Created by rja on 06/02/2021.
//

#include "FieldX.h"
#include "src/core/table/TableBaseX.h"

FieldBaseX::FieldBaseX(RowX *row, bool is_key, size_t nelement, const std::type_info &type_info, size_t element_size) :
        m_row(row), m_is_key(is_key), m_element_size(element_size), m_nelement(nelement),
        m_size(nelement*element_size), m_type_info(type_info), m_offset(row->add_field(this)){
    m_details["number of elements"] = std::to_string(m_nelement);
    m_details["field size (bytes)"] = std::to_string(m_size);
    m_details["offset from row (bytes)"] = std::to_string(m_offset);
}

FieldBaseX::FieldBaseX(const FieldBaseX &other) :
        FieldBaseX(other.m_row->m_last_copied, other.m_is_key, other.m_nelement, other.m_type_info, other.m_element_size){}

char *FieldBaseX::begin() {
    return (char*)(m_row->m_dbegin)+m_offset;
}

const char *FieldBaseX::begin() const {
    return (const char*)(m_row->m_dbegin)+m_offset;
}

char *FieldBaseX::raw_view() {
    ASSERT(m_row->m_dbegin);
    ASSERT(m_element_offset<m_size);
    return (char*)(m_row->m_dbegin)+m_offset+m_element_offset;
}

const char *FieldBaseX::raw_view() const {
    ASSERT(m_row->m_dbegin);
    ASSERT(m_element_offset<m_size);
    return (const char*)(m_row->m_dbegin)+m_offset+m_element_offset;
}

void FieldBaseX::zero() {
    std::memset(raw_view(), 0, m_element_size);
}

void FieldBaseX::zero_all() {
    std::memset((m_row->m_dbegin)+m_offset, 0, m_size);
}

RowX::RowX(const RowX &other) : m_table(other.m_table),
    m_dbegin(m_table?m_table->dbegin(): nullptr) {
    other.m_last_copied = this;
}

bool RowX::select_first() const {
    ASSERT(m_table);
    m_dbegin = m_table->dbegin();
    return true;
}

bool RowX::select_next() const {
    ASSERT(m_table);
    m_dbegin += m_dsize;
    if(m_dbegin>=m_table->m_bw.m_dend) {
        // Row accesses Table out of bounds
        select_first();
        return false;
    }
    return true;
}

bool RowX::select(const size_t &irow) const {
    ASSERT(m_table);
    m_dbegin = m_table->dbegin()+irow*m_dsize;
    if(m_dbegin>=m_table->m_bw.m_dend) {
        // Row accesses Table out of bounds
        select_first();
        return false;
    }
    return true;
}
