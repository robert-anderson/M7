//
// Created by rja on 09/02/2021.
//

#include "FieldBase.h"
#include "src/core/parallel/MPIAssert.h"


FieldBase::FieldBase(Row* row, size_t item_size, size_t nitem, const std::type_info &type_info) :
        m_item_size(item_size), m_type_info(type_info), m_nitem(nitem),
        m_size(m_item_size*m_nitem), m_null_string(m_size, 0){
    m_row = row;
    if (m_row) m_row_offset = m_row->add_field(this);
}

FieldBase::FieldBase(const FieldBase &other) :
        FieldBase(other.m_row ? other.m_row->m_child : nullptr, other.m_item_size, other.m_nitem, other.m_type_info) {}

bool FieldBase::is_comparable(const FieldBase &other) const {
    return m_item_size==other.m_item_size &&
           m_type_info==other.m_type_info && m_size==other.m_size;
}

FieldBase &FieldBase::operator=(const FieldBase &other) {
    if (&other == this) return *this;
    MPI_ASSERT(is_comparable(other),
               "can't copy from a field which is either incompatible or has a different selection length")
    std::memcpy(begin(), other.begin(), m_size);
    return *this;
}

void FieldBase::add_to_row(Row *row) {
    if (!row) return;
    MPI_ASSERT(m_size, "Can't add a field of zero size to a row.");
    MPI_REQUIRE(!is_added_to_row(), "Field must not be already associated with a row");
    m_row_offset = row->add_field(this);
    m_row = row;
}

bool FieldBase::is_added_to_row() const {
    return (m_row_offset != ~0ul) && m_row;
}

char *FieldBase::begin() const {
    MPI_ASSERT(is_added_to_row(), "Field is not associated with row");
    MPI_ASSERT(m_row->m_dbegin, "Associated row is not pointing to table buffer");
    return (char *) (m_row->m_dbegin) + m_row_offset;
}

char *FieldBase::begin(const size_t &iitem) const {
    return begin()+iitem*m_item_size;
}

char *FieldBase::end() const {
    return begin()+m_size;
}
char *FieldBase::end(const size_t &iitem) const {
    return end() + iitem * m_item_size;
}

Row *FieldBase::row_of_copy() const {
    return m_row? m_row->m_child : nullptr;
}

void FieldBase::zero() {
    std::memset(begin(), 0, m_size);
}

bool FieldBase::is_zero() const {
    ASSERT(!m_null_string.empty());
    return std::memcmp(begin(), m_null_string.data(), m_size) == 0;
}

bool FieldBase::operator==(const FieldBase &other) const {
    MPI_ASSERT(is_comparable(other),
               "can't copy from a field which is either incompatible or has a different selection length")
    return std::memcmp(begin(), other.begin(), m_size) == 0;
}

defs::hash_t FieldBase::hash() const {
    return hashing::fnv_hash(begin(), m_size);
}

std::string FieldBase::to_string() const {
    std::string tmp;
    for (size_t i=0ul; i<m_nitem; ++i) tmp+=to_string_element(i);
    return tmp;
}