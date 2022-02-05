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

StringField::StringField(Row *row, size_t length, std::string name) : base_t(row, {length}, name) {}

StringField::StringField(const StringField &other) : base_t(other){}

StringField &StringField::operator=(const StringField &other) {
    base_t::operator=(other);
    return *this;
}

StringField &StringField::operator=(const char *str) {
    auto len = std::strlen(str);
    DEBUG_ASSERT_LE(len, nelement(), "String length does not match that of string field");
    zero();
    std::copy(str, str+len, dbegin());
    return *this;
}

std::string StringField::to_string() const {
    auto len = std::min(nelement(), std::strlen(dbegin()));
    return {dbegin(), len};
}
