//
// Created by Robert J. Anderson on 09/02/2021.
//

#include "Row.h"
#include "FieldBase.h"


str_t Row::field_names_string() const {
    str_t tmp;
    uint_t i = 0ul;
    for (const FieldBase* field : m_fields)
        tmp+="field "+std::to_string(i++)+": " + field->m_name + "\n";
    return tmp;
}

str_t Row::to_string() const {
    str_t tmp;
    for (auto field: m_fields) tmp += " " + field->to_string();
    return tmp;
}

uint_t Row::add_field(FieldBase *field) {
    REQUIRE_TRUE(field, "Field pointer should not be null");
    // returns the offset in bytes for the column being added
    auto offset = 0ul;
    if (!m_fields.empty()) {
        offset = m_fields.back()->m_row_offset + m_fields.back()->m_size;
        if (!(m_fields.back()->m_type_info == field->m_type_info)) {
            // go to next whole system word
            offset = integer::divceil(offset, Buffer::c_nbyte_word) * Buffer::c_nbyte_word;
        }
    }

    m_current_offset = offset + field->m_size;
    m_size = integer::divceil(m_current_offset, Buffer::c_nbyte_word) * Buffer::c_nbyte_word;

    m_fields.push_back(field);
    return offset;
}

uint_t Row::nfield() const {
    return m_fields.size();
}

void Row::clear() {
    std::fill(m_begin, m_begin+m_size, 0);
}

bool Row::is_cleared() const {
    return m_table->is_cleared(index());
}

bool Row::is_h5_write_exempt() const {
    return false;
}

bool Row::is_protected() const {
    return m_table->is_protected(m_i);
}

Row::Row(const Row &other) {
    m_table = other.m_table;
    m_begin = other.m_begin;
    other.m_child = this;
    REQUIRE_TRUE(m_fields.empty(), "newly copied Row has associated Fields");
}

Row::~Row() {}
