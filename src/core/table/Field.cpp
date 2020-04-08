//
// Created by Robert John Anderson on 2020-03-26.
//

#include <climits>
#include "Field.h"
#include "Table.h"
#include "Element.h"

Field::Field(Table *table, size_t element_size, size_t nelement, const std::type_info &type_info,
    const std::string& description):
m_table(table), m_element_size(element_size),
m_element_dsize(integer_utils::divceil(element_size, sizeof(defs::data_t))),
m_nelement(nelement), m_type_index(type_info), m_description(description),
m_offset(table ? table->add_field(this) : ~0ul){}

char *Field::begin(const size_t &irow, const size_t &isegment) {
    return m_table->field_begin(this, irow, isegment);
}

Element Field::operator()(const size_t &irow, const size_t &isegment, const size_t &ielement) {
    return Element(this, element_begin(irow, isegment, ielement));
}

size_t Field::nbit() const{
    return CHAR_BIT*m_element_size;
}

size_t Field::element_dsize() const {
    return m_element_dsize;
}

bool Field::compatible_with(const Field &rhs) const {
    return m_element_size == rhs.m_element_size &&
           m_nelement == rhs.m_nelement &&
           m_type_index == rhs.m_type_index;
}

char *Field::element_begin(const size_t &irow, const size_t &isegment, const size_t &ielement) {
    return begin(irow, isegment) + ielement * m_element_size;
}

void Field::expand_table(size_t delta_nrow) {
    m_table->expand(delta_nrow);
}

bool Field::is_allocated() const {
    return m_table->is_allocated();
}
