//
// Created by Robert J. Anderson on 09/02/2021.
//

#include <M7_lib/parallel/MPIAssert.h>

#include <utility>
#include "FieldBase.h"


FieldBase::FieldBase(Row *row, uint_t size, const std::type_info &type_info, str_t name, bool force_own_words) :
        m_type_info(type_info), m_size(size), m_name(std::move(name)),
        m_force_own_words(force_own_words), m_null_string(std::max(1ul, m_size), 0) {
    if (!row) return;
    REQUIRE_FALSE(belongs_to_row(), "Field must not be already associated with a row");
    m_row_offset = row->add_field(this);
    m_row = row;
    m_row_index = m_row->m_fields.size()-1;
}

FieldBase::FieldBase(const FieldBase &other) :
        FieldBase(other.row_of_copy(), other.m_size, other.m_type_info, other.m_name, other.m_force_own_words) {}

bool FieldBase::is_comparable(const FieldBase &other) const {
    return m_type_info == other.m_type_info && m_size == other.m_size;
}

bool FieldBase::belongs_to_row() const {
    return (m_row_offset != ~0ul) && (m_row_index != ~0ul) && m_row;
}

bool FieldBase::belongs_to_row(const Row* row) const {
    return m_row==row;
}

bool FieldBase::belongs_to_row(const Row& row) const {
    return belongs_to_row(&row);
}

buf_t *FieldBase::begin() const {
    DEBUG_ASSERT_TRUE(belongs_to_row(), "Field is not associated with row");
    return m_row->begin() + m_row_offset;
}

buf_t *FieldBase::end() const {
    return begin() + m_size;
}

const Row *FieldBase::row() const {
    return m_row;
}

Row *FieldBase::row_of_copy() const {
    if (!m_row) return nullptr;
    REQUIRE_TRUE(m_row->m_child, "Row did not produce a copy");
    return m_row->m_child;
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

hash::digest_t FieldBase::hash() const {
    return hash::fnv(begin(), m_size);
}

uint_t FieldBase::to_buffer(buf_t* buf, uint_t irow_begin, uint_t nitem_max, std::set<uint_t> irows_empty) const {
    DEBUG_ASSERT_LT(irow_begin, m_row->m_table->nrow_in_use(), "row index OOB");
    uint_t nitem = 0ul;
    for (uint_t irow = irow_begin; irow<m_row->m_table->nrow_in_use(); ++irow){
        if (nitem==nitem_max) return nitem;
        if (irows_empty.find(irow)!=irows_empty.end()) continue;
        auto dst = buf + irow*m_size;
        auto src = m_row->m_table->begin()+m_row->m_size*irow+m_row_offset;
        std::memcpy(dst, src, m_size);
    }
    return nitem;
}
