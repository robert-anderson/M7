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
        const auto prev = m_fields.back();
        offset = prev->m_row_offset + prev->m_size;
        if (field->m_force_own_words || prev->m_force_own_words || (prev->m_type_info != field->m_type_info)) {
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

void Row::free() {
    m_table->free(index());
}

bool Row::is_h5_write_exempt() const {
    return false;
}

Row::Row(const Row &other) {
    m_table = other.m_table;
    m_begin = other.m_begin;
    other.m_child = this;
    REQUIRE_TRUE(m_fields.empty(), "newly copied Row has associated Fields");
}

Row::~Row() {}

void Row::protect() {
    m_table->protect(index());
}

uint_t Row::protection_level() const {
    return m_table->protection_level(index());
}

bool Row::is_protected() const {
    return m_table->is_protected(index());
}

void Row::release() {
    m_table->release(index());
}
