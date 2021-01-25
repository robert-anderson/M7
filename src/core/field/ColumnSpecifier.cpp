//
// Created by rja on 21/10/2020.
//

#include "ColumnSpecifier.h"
#include "src/core/table/Table.h"

/*
 * View
 */
ColumnSpecifier::View::View(const ColumnSpecifier &field, char *ptr) : m_spec(field), m_ptr(ptr){}

ColumnSpecifier::View::View(const ColumnSpecifier::View &other) : m_spec(other.m_spec), m_ptr(other.m_ptr){}

const size_t &ColumnSpecifier::View::element_size() const {
    return m_spec.element_size();
}

int ColumnSpecifier::View::compare(const ColumnSpecifier::View &other) const {
    ASSERT(element_size()==other.element_size());
    return std::memcmp(m_ptr, other.m_ptr, element_size());
}

bool ColumnSpecifier::View::operator==(const ColumnSpecifier::View &other) const {
    return compare(other)==0;
}

bool ColumnSpecifier::View::operator!=(const ColumnSpecifier::View &other) {
    return !(*this==other);
}

bool ColumnSpecifier::View::is_zero() const {
    return !std::memcmp(m_ptr, m_spec.m_null_buffer.data(), element_size());
}

void ColumnSpecifier::View::print() const {
    std::cout << to_string() << std::endl;
}

const defs::data_t *ColumnSpecifier::View::cdptr(const size_t &i) const {
    ASSERT(i * defs::nbyte_data < m_spec.element_size());
    return ((defs::data_t *) m_ptr) + i;
}

defs::data_t *ColumnSpecifier::View::dptr(const size_t &i) const {
    ASSERT(i * defs::nbyte_data < m_spec.element_size());
    return ((defs::data_t *) m_ptr) + i;
}

void ColumnSpecifier::View::zero() {
    std::memset(m_ptr, 0, element_size());
}

void ColumnSpecifier::View::mpi_bcast(size_t iroot) {
    mpi::bcast(m_ptr, element_size(), iroot);
}

ColumnSpecifier::View &ColumnSpecifier::View::operator=(const ColumnSpecifier::View &other) {
    ASSERT(m_spec.element_size() == other.m_spec.element_size());
    if (&other != this)
        std::memcpy(m_ptr, other.m_ptr, m_spec.element_size());
    return *this;
}

/*
 * ColumnSpecifier
 */
ColumnSpecifier::ColumnSpecifier(size_t element_size, const std::type_info &type_info) :
        m_data{element_size, type_info, {{"element size (bytes)", std::to_string(element_size)}}},
        m_null_buffer(element_size, 0){}

const size_t &ColumnSpecifier::element_size() const {
    return m_data.m_element_size;
}

const std::type_info &ColumnSpecifier::type_info() const {
    return m_data.m_type_info;
}

defs::hash_t ColumnSpecifier::hash(const ColumnSpecifier::View &view) {
    return hashing::fnv_hash(view.m_ptr, view.element_size());
}

bool ColumnSpecifier::comparable_with(const ColumnSpecifier &other) const {
    return (type_info()==other.type_info()) && (element_size()==other.element_size());
}