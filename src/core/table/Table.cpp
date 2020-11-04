//
// Created by rja on 21/10/2020.
//

#include "Table.h"

size_t TableX::push_back() {
    if (m_hwm>=m_nrow) throw std::runtime_error("Table capacity reached");
    return m_hwm++;
}

char *TableX::begin() {
    return (char *) m_bw.m_ptr;
}

char *TableX::begin(const size_t &irow) {
    ASSERT(irow<m_hwm)
    return begin() + irow * m_row_size;
}

size_t TableX::add_field(const TableField *field) {
    // returns the offset in bytes for the field being added
    auto offset = 0ul;
    if(!m_fields.empty()){
        offset = m_fields.back()->m_offset+m_fields.back()->m_size;
        if (!m_fields.back()->is_same_type_as(*field)){
            // go to next whole dataword
            offset = integer_utils::divceil(offset, defs::nbyte_data)*defs::nbyte_data;
        }
    }

    m_tight_row_size = offset+field->m_size;
    m_row_dsize = integer_utils::divceil(m_tight_row_size, defs::nbyte_data);
    m_row_size = m_row_dsize*defs::nbyte_data;

    m_fields.push_back(field);
    return offset;
}

void TableX::move(BufferWindow new_bw) {
    if (m_bw) std::memmove(m_bw.m_ptr, new_bw.m_ptr, sizeof(defs::data_t) * std::min(m_bw.m_dsize, new_bw.m_dsize));
    m_bw = new_bw;
    if (!m_row_size) return;
    m_nrow = (sizeof(defs::data_t) * m_bw.m_dsize) / m_row_size;
}

void TableX::clear() {
    std::memset((char *) (m_bw.m_ptr), 0, m_row_size * m_hwm);
    m_hwm = 0ul;
}

void TableX::clear_row(const size_t &irow) {
    std::memset(begin(irow), 0, m_tight_row_size);
}

std::string TableX::field_details(size_t width) const {
    std::string res;
    for (size_t i = 0ul; i < m_fields.size(); ++i) {
        std::string desc;
        if (!m_fields[i]->m_description.empty()) desc = " (\""+m_fields[i]->m_description+"\")";
        res += "\nField " + std::to_string(i) + desc + ":\n";
        for (auto pair:m_fields[i]->m_data.m_details)
            res += "\t" + utils::padded_string(pair.first, width) + ": " + utils::padded_string(pair.second, width)+"\n";
    }
    return res;
}

void TableX::print_field_details(size_t width) const {
    std::cout << field_details(width) << std::endl;
}
