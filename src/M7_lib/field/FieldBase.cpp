//
// Created by Robert J. Anderson on 09/02/2021.
//

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/hdf5/Field.h>

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
    DEBUG_ASSERT_EQ(m_null_string.size(), m_size, "comparison-with-zero string is incorrect length");
    return std::memcmp(cbegin(), m_null_string.data(), m_size) == 0;
}

int FieldBase::cmp(const FieldBase &other) const {
    return std::memcmp(cbegin(), other.cbegin(), m_size);
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
    return hash::fnv(cbegin(), m_size);
}

void FieldBase::save_fn(const hdf5::NodeWriter& nw, const str_t& name, bool this_rank, uint_t max_nitem_per_op) const {
    hdf5::field::save<buf_t>(*this, nw, name, {m_size}, {"bytes"}, max_nitem_per_op, this_rank);
}

void FieldBase::load_fn(const hdf5::NodeReader& nr, const str_t& name,
                        uint_t max_nitem_per_op, bool part, bool this_rank) {
    hdf5::field::load<buf_t>(*this, nr, name, max_nitem_per_op, part, this_rank);
}

void FieldBase::save(const hdf5::NodeWriter& nw, const str_t& name, bool this_rank, uint_t max_nitem_per_op) const {
    save_fn(nw, name, this_rank, max_nitem_per_op);
}

void FieldBase::save(const hdf5::NodeWriter& nw, const str_t& name, bool this_rank) const {
    save(nw, name, hdf5::field::c_default_max_nitem_per_op, this_rank);
}

void FieldBase::save(const hdf5::NodeWriter& nw, bool this_rank) const {
    save(nw, m_name, hdf5::field::c_default_max_nitem_per_op, this_rank);
}

void FieldBase::load(const hdf5::NodeReader& nr, const str_t& name, bool part, bool this_rank, uint_t max_nitem_per_op) {
    load_fn(nr, name, max_nitem_per_op, part, this_rank);
}

void FieldBase::load(const hdf5::NodeReader& nr, const str_t& name, bool part, bool this_rank) {
    load(nr, name, hdf5::field::c_default_max_nitem_per_op, part, this_rank);
}

void FieldBase::load(const hdf5::NodeReader& nr, bool part, bool this_rank) {
    load(nr, m_name, part, this_rank);
}
