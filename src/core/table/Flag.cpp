//
// Created by Robert John Anderson on 2020-03-26.
//

#include "Flag.h"
#include "FlagField.h"


Flag::Flag(FlagField* field, size_t nelement, const std::string &description):
m_field(field), m_nelement(nelement), m_description(description), m_offset(field->add_flag(this)) {}

FlagElement Flag::operator()(const size_t &irow, const size_t &isegment) {
    return FlagElement(this, (*m_field)(irow, isegment));
}

FlagElement::FlagElement(Flag *flag, BitsetElement bitset_element) :
    m_flag(flag), m_bitset_element(bitset_element){}

void FlagElement::operator=(bool v) {
    if (v){
        m_bitset_element.set(m_flag->offset());
    } else {
        m_bitset_element.clr(m_flag->offset());
    }
}

FlagElement::operator bool() {
    return m_bitset_element.get(m_flag->offset());
}
