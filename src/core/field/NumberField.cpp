//
// Created by anderson on 1/26/22.
//

#include "NumberField.h"

NumberFieldBase::NumberFieldBase(Row *row, size_t element_size, size_t nelement, bool is_complex,
                                 const std::type_info &type_info, std::string name) :
        FieldBase(row, element_size * nelement, type_info, name),
        m_element_size(element_size), m_nelement(nelement), m_is_complex(is_complex) {}

NumberFieldBase::NumberFieldBase(const NumberFieldBase &other) :
    FieldBase(other), m_element_size(other.m_element_size),
    m_nelement(other.m_nelement), m_is_complex(other.m_is_complex){}

NumberFieldBase::NumberFieldBase(NumberFieldBase &&other) :
    FieldBase(std::move(other)),  m_element_size(other.m_element_size),
    m_nelement(other.m_nelement), m_is_complex(other.m_is_complex){}

NumberFieldBase &NumberFieldBase::operator=(NumberFieldBase &&other) {
    *this = other;
    return *this;
}
