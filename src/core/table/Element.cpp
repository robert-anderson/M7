//
// Created by Robert John Anderson on 2020-03-26.
//

#include <src/consts.h>
#include <src/core/hash/Hashing.h>
#include <src/utils.h>
#include "Element.h"
#include "Field.h"

Element::Element(Field *field, char *begin) : m_field(field), m_begin(begin){}

size_t Element::hash() const{
    return hashing::fnv_hash(m_begin, size());
}

size_t Element::size() const {
    return m_field->m_element_size;
}

size_t Element::dsize() const {
    return m_field->m_element_dsize;
}

std::string Element::to_string() const {
    return std::string(m_begin, size());
}

bool Element::compatible_with(const Element &rhs) const {
    return m_field->compatible_with(*rhs.m_field);
}

defs::data_t &Element::dataword(const size_t &idataword) const {
    return *(((defs::data_t*)m_begin)+idataword);
}

defs::data_t Element::get_dataword(const size_t &idataword) const {
    return dataword(idataword);
}

defs::data_t Element::get_dataword(const size_t &idataword, const size_t &nbit) const {
    defs::data_t tmp = get_dataword(idataword);
    return bit_utils::truncate(tmp, nbit);
}

defs::data_t Element::get_antidataword(const size_t &idataword) const {
    return ~*(((defs::data_t*)m_begin)+idataword);
}

defs::data_t Element::get_antidataword(const size_t &idataword, const size_t &nbit) const {
    defs::data_t tmp = get_antidataword(idataword);
    return bit_utils::truncate(tmp, nbit);
}

int Element::cmp(const Element &rhs) const {
    return std::memcmp(m_begin, rhs.m_begin, m_field->m_element_size);
}

Element& Element::operator=(const Element& rhs) {
    assert(compatible_with(rhs));
    std::memcpy(m_begin, rhs.m_begin, size());
    return *this;
}

bool Element::operator==(const Element &rhs) const {
    return cmp(rhs)==0;
}

bool Element::operator!=(const Element &rhs) const {
    return cmp(rhs)!=0;
}

char *Element::begin() const {
    return m_begin;
}

size_t Element::nbit() const {
    return m_field->nbit();
}

void Element::zero() {
    std::memset(m_begin, 0, size());
}

bool Element::is_zero() const {
    std::vector<char> zero(m_field->m_element_size, 0);
    return memcmp(m_begin, zero.data(), m_field->m_element_size);
}

void Element::print() const {
    std::cout << to_string() << std::endl;
}

const Field *Element::field() const {return m_field;}