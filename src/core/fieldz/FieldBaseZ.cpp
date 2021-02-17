//
// Created by rja on 09/02/2021.
//

#include "FieldBaseZ.h"
#include "src/core/parallel/MPIAssert.h"


FieldBaseZ::FieldBaseZ(RowZ* row, size_t item_size, size_t nitem, const std::type_info &type_info) :
        m_item_size(item_size), m_type_info(type_info), m_nitem(nitem), m_size(m_item_size*m_nitem){
    m_row = row;
    if (m_row) m_row_offset = m_row->add_field(this);
}

FieldBaseZ::FieldBaseZ(const FieldBaseZ &other) :
    FieldBaseZ(other.m_row, other.m_item_size, other.m_nitem, other.m_type_info) {
    m_row_offset = other.m_row_offset;
}

bool FieldBaseZ::is_added_to_row() const {
    return (m_nitem != ~0ul) && (m_size != ~0ul) && (m_row_offset != ~0ul) && m_row;
}

char *FieldBaseZ::begin() const {
    MPI_ASSERT(is_added_to_row(), "Field is not associated with row");
    MPI_ASSERT(m_row->m_dbegin, "Associated row is not pointing to table buffer");
    return (char *) (m_row->m_dbegin) + m_row_offset;
}

char *FieldBaseZ::begin(const size_t &iitem) const {
    return begin()+iitem*m_item_size;
}

char *FieldBaseZ::end() const {
    return begin()+m_size;
}
char *FieldBaseZ::end(const size_t &iitem) const {
    return end()+iitem*m_item_size;
}

void FieldBaseZ::zero() {
    std::memset(begin(), 0, m_size);
}

bool FieldBaseZ::is_zero() const {
    return std::memcmp(begin(), m_null_item_string.data(), m_size) == 0;
}

defs::hash_t FieldBaseZ::hash() const {
    return hashing::fnv_hash(begin(), m_size);
}

std::string FieldBaseZ::to_string() const {
    std::string tmp;
    for (size_t i=0ul; i<m_nitem; ++i) tmp+=to_string_element(i);
    return tmp;
}