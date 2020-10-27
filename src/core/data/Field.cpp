//
// Created by rja on 21/10/2020.
//

#include "Field.h"
#include "Table.h"

FieldBaseX::FieldBaseX(size_t element_size, const std::type_info &type_info) :
        m_element_size(element_size), m_type_info(type_info) {}

bool FieldBaseX::is_same_type_as(const FieldBaseX &other) const {
    return m_type_info == other.m_type_info;
}
