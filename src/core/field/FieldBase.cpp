//
// Created by rja on 09/02/2021.
//

#include "FieldBase.h"
#include "src/core/parallel/MPIAssert.h"


FieldBase::FieldBase(Row *row, size_t size, const std::type_info &type_info, std::string name) :
        m_type_info(type_info), m_size(size),
        m_name(name), m_null_string(std::max(1ul, m_size), 0) {
    m_row = row;
    if (m_row) m_row_offset = m_row->add_field(this);
}

FieldBase::FieldBase(const FieldBase &other) :
        FieldBase(row_of_copy(), other.m_size, other.m_type_info, other.m_name) {}

bool FieldBase::is_comparable(const FieldBase &other) const {
    return m_type_info == other.m_type_info && m_size == other.m_size;
}

void FieldBase::add_to_row(Row *row) {
    if (!row) return;
    REQUIRE_FALSE(belongs_to_row(), "Field must not be already associated with a row");
    m_row_offset = row->add_field(this);
    m_row = row;
}

bool FieldBase::belongs_to_row() const {
    return (m_row_offset != ~0ul) && m_row;
}

bool FieldBase::belongs_to_row(const Row* row) const {
    return m_row==row;
}

bool FieldBase::belongs_to_row(const Row& row) const {
    return belongs_to_row(&row);
}

char *FieldBase::begin() const {
    DEBUG_ASSERT_TRUE(belongs_to_row(), "Field is not associated with row");
    return m_row->begin() + m_row_offset;
}

char *FieldBase::end() const {
    return begin() + m_size;
}

Row *FieldBase::row_of_copy() const {
    return m_row ? m_row->m_child : nullptr;
}

void FieldBase::zero() {
    std::memset(begin(), 0, m_size);
}

bool FieldBase::is_zero() const {
    ASSERT(!m_null_string.empty());
    return std::memcmp(begin(), m_null_string.data(), m_size) == 0;
}

int FieldBase::cmp(const FieldBase &other) const {
    return std::memcmp(begin(), other.begin(), m_size);
}

bool FieldBase::operator==(const FieldBase &other) const {
    return cmp(other) == 0;
}

bool FieldBase::operator!=(const FieldBase &other) const {
    return cmp(other) != 0;
}

bool FieldBase::operator<(const FieldBase &other) const {
    return cmp(other) < 0;
}

bool FieldBase::operator>(const FieldBase &other) const {
    return cmp(other) > 0;
}

bool FieldBase::operator<=(const FieldBase &other) const {
    return cmp(other) <= 0;
}

bool FieldBase::operator>=(const FieldBase &other) const {
    return cmp(other) >= 0;
}

defs::hash_t FieldBase::hash() const {
    return hashing::fnv_hash(begin(), m_size);
}

FieldBase &FieldBase::operator=(const FieldBase &other) {
    if (&other == this) return *this;
    DEBUG_ASSERT_TRUE(is_comparable(other), "can't compare to incompatible field");
    std::memcpy(begin(), other.begin(), m_size);
    return *this;
}