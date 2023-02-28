//
// Created by Robert J. Anderson on 1/26/22.
//

#include "NumberField.h"

#include <utility>

NumberFieldBase::NumberFieldBase(Row *row, uint_t element_size, uint_t nelement, bool is_complex,
                                 const std::type_info &type_info, str_t name, bool force_own_words) :
        FieldBase(row, element_size * nelement, type_info, std::move(name), force_own_words),
        m_element_size(element_size), m_nelement(nelement), m_is_complex(is_complex) {}

NumberFieldBase::NumberFieldBase(const NumberFieldBase &other) :
    FieldBase(other), m_element_size(other.m_element_size),
    m_nelement(other.m_nelement), m_is_complex(other.m_is_complex){}

StringField::StringField(Row *row, uint_t length, const str_t& name, bool force_own_words) :
    base_t(row, {length}, name, force_own_words) {}

StringField::StringField(const StringField &other) : base_t(other){}

StringField &StringField::operator=(const StringField &other) {
    base_t::operator=(other);
    return *this;
}

StringField &StringField::operator=(const char *str) {
    auto len = std::strlen(str);
    DEBUG_ASSERT_LE(len, nelement(), "String length does not match that of string field");
    zero();
    std::copy(str, str+len, tbegin());
    return *this;
}

StringField &StringField::operator=(const str_t &str) {
    *this=str.c_str();
    return *this;
}

bool StringField::operator==(const char *str) const {
    return !memcmp(str, ctbegin(), std::strlen(str));
}

bool StringField::operator==(const str_t &str) const {
    return (*this==str.c_str());
}

bool StringField::operator!=(const char *str) const {
    return !(*this==str);
}

bool StringField::operator!=(const str_t &str) const {
    return !(*this==str);
}

str_t StringField::to_string() const {
    auto len = std::min(nelement(), std::strlen(ctbegin()));
    return {ctbegin(), len};
}
