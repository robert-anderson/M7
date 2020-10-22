//
// Created by rja on 02/10/2020.
//

#include "FieldBase.h"

FieldBase::FieldBase(Table_NEW *table, size_t element_size, size_t nelement, const std::type_info &type_info, std::string description):
        m_table(table), m_element_size(element_size),
        m_nelement(nelement), m_size(m_element_size*m_nelement),
        m_dsize(integer_utils::divceil(m_size, sizeof(defs::data_t))),
        m_type_info(type_info), m_description(std::move(description)){}

bool FieldBase::is_same_type_as(const FieldBase &other) const {
    return other.m_type_info == m_type_info;
}

char *FieldBase::begin(const size_t &irow) const {
    return m_table->row_ptr(irow)+m_offset;
}

char *FieldBase::begin(const size_t &irow, const size_t &ielement) const {
    return m_table->row_ptr(irow)+m_offset+ielement*m_element_size;
}

void FieldBase::set_offsets() {
    if (!m_table->m_fields.empty()) {
        const auto &last_field = *m_table->m_fields.back();
        if (is_same_type_as(last_field)) m_offset = last_field.back_offset();
        else m_offset = integer_utils::round_up(last_field.back_offset(), sizeof(defs::data_t));
    }
    m_table->m_tight_row_size = back_offset();
    m_table->add_field(this);
}
