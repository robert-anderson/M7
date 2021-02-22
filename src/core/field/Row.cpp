//
// Created by rja on 09/02/2021.
//

#include "Row.h"
#include "FieldBase.h"

std::string Row::to_string() const {
    std::string tmp;
    for (auto field: m_fields) tmp += " " + field->to_string();
    return tmp;
}

size_t Row::add_field(FieldBase *field) {
    // returns the offset in bytes for the column being added
    auto offset = 0ul;
    if (!m_fields.empty()) {
        offset = m_fields.back()->m_row_offset + m_fields.back()->m_size;
        if (!(m_fields.back()->m_type_info == field->m_type_info)) {
            // go to next whole dataword
            offset = integer_utils::divceil(offset, defs::nbyte_data) * defs::nbyte_data;
        }
    }

    m_current_offset = offset + field->m_size;
    m_dsize = integer_utils::divceil(m_current_offset, defs::nbyte_data);
    m_size = m_dsize * defs::nbyte_data;

    m_fields.push_back(field);
    return offset;
}

void Row::clear() {
    std::fill(m_dbegin, m_dbegin+m_dsize, 0);
}

bool Row::is_cleared() const {
    return std::all_of(m_dbegin, m_dbegin+m_dsize, [](const defs::data_t& d){return d==0;});
}
