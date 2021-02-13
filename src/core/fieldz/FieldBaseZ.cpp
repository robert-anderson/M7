//
// Created by rja on 09/02/2021.
//

#include "FieldBaseZ.h"
#include "src/core/parallel/MPIAssert.h"

FieldBaseZ::FieldBaseZ(size_t element_size, const std::type_info &type_info) :
        m_element_size(element_size), m_type_info(type_info) {}

FieldBaseZ::FieldBaseZ(const FieldBaseZ &other) : FieldBaseZ(other.m_element_size, other.m_type_info) {}

bool FieldBaseZ::is_added_to_row() const {
    return (m_nelement != ~0ul) && (m_size != ~0ul) && (m_row_offset != ~0ul) && m_row;
}

void FieldBaseZ::restart() const {
    MPI_ASSERT(is_added_to_row(), "Field isn't added to a row");
    m_view_offset = m_row_offset;
    MPI_ASSERT(!oob(), "Field has jumped to out of bounds view")
}

void FieldBaseZ::step() const {
    ASSERT(is_added_to_row());
    m_view_offset += m_element_size;
    MPI_ASSERT(!oob(), "Field has jumped to out of bounds view")
}

void FieldBaseZ::jump(const size_t &iflat) const {
    ASSERT(is_added_to_row());
    m_view_offset = m_row_offset + m_element_size * iflat;
    MPI_ASSERT(!oob(), "Field has jumped to out of bounds view")
}

bool FieldBaseZ::try_restart() const {
    ASSERT(is_added_to_row());
    m_view_offset = m_row_offset;
    return !oob();
}

bool FieldBaseZ::try_step() const {
    ASSERT(is_added_to_row());
    m_view_offset += m_element_size;
    if (oob()) {
        m_view_offset = m_row_offset;
        return false;
    }
    return true;
}

bool FieldBaseZ::try_jump(const size_t& iflat) const {
    ASSERT(is_added_to_row());
    m_view_offset = m_row_offset + m_element_size * iflat;
    if (oob()) {
        m_view_offset = m_row_offset;
        return false;
    }
    return true;
}

FieldBaseZ &FieldBaseZ::copy(const FieldBaseZ &other) {
    MPI_REQUIRE(m_element_size==other.m_element_size &&
                m_type_info==other.m_type_info, "incompatible formats");
    std::memcpy(raw_view(), other.raw_view(), m_element_size);
    return *this;
}

FieldBaseZ &FieldBaseZ::copy_all(const FieldBaseZ &other) {
    std::memcpy(begin(), other.begin(), m_size);
    return *this;
}

//FieldBaseZ &FieldBaseZ::operator=(const FieldBaseZ &other) {
//    copy(other);
//    return *this;
//}

char *FieldBaseZ::begin() const {
    MPI_ASSERT(is_added_to_row(), "Field is not associated with row");
    MPI_ASSERT(m_row->m_dbegin, "Associated row is not pointing to table buffer");
    return (char *) (m_row->m_dbegin) + m_row_offset;
}

char *FieldBaseZ::raw_view() const {
    MPI_ASSERT(is_added_to_row(), "Field is not associated with row");
    MPI_ASSERT(m_row->m_dbegin, "Associated row is not pointing to table buffer");
    return (char *) (m_row->m_dbegin) + m_view_offset;
}

void FieldBaseZ::zero() {
    std::memset(raw_view(), 0, m_element_size);
}

void FieldBaseZ::zero_all() {
    std::memset(begin(), 0, m_size);
}

bool FieldBaseZ::is_zero() const {
    return std::memcmp(raw_view(), m_null_element_string.data(), m_element_size) == 0;
}

bool FieldBaseZ::equals(const FieldBaseZ &other) const {
    return std::memcmp(raw_view(), other.raw_view(), m_size) == 0;
}

bool FieldBaseZ::is_same_type_as(const FieldBaseZ &other) const {
    return m_type_info == other.m_type_info;
}

defs::hash_t FieldBaseZ::hash() const {
    return hashing::fnv_hash(begin(), m_size);
}

std::string FieldBaseZ::to_string_all() const {
    std::string tmp;
    bool loop = try_restart();
    while (loop) {
        tmp += to_string() + " ";
        loop = try_step();
    }
    return tmp;
}