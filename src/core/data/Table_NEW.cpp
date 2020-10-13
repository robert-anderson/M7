//
// Created by rja on 02/10/2020.
//

#include <src/core/util/utils.h>
#include "Table_NEW.h"
#include "FieldBase.h"

Table_NEW::Table_NEW(BufferWindow bw) : m_bw(std::move(bw)){}

Table_NEW::Table_NEW() : Table_NEW(BufferWindow()){}

void Table_NEW::move(BufferWindow new_bw) {
    if (m_bw) std::memmove(m_bw.m_ptr, new_bw.m_ptr, sizeof(defs::data_t) * std::min(m_bw.m_dsize, new_bw.m_dsize));
    m_bw = new_bw;
    m_nrow = (sizeof(defs::data_t)*m_bw.m_dsize)/m_row_size;
}

size_t Table_NEW::push_back() {
    ASSERT(m_hwm<m_nrow);
    return m_hwm++;
}

void Table_NEW::add_field(FieldBase *field) {
    m_ncacheline = integer_utils::divceil(m_tight_row_size, defs::ncacheline_byte);
    m_row_size = m_ncacheline*defs::ncacheline_byte;
    m_fields.push_back(field);
}

char *Table_NEW::row_ptr(const size_t &irow) const {
    ASSERT(irow<m_hwm);
    ASSERT(m_bw.m_dsize);
    return (char*)(m_bw.m_ptr)+irow*m_row_size;
}

void Table_NEW::clear() {
    std::memset((char*)(m_bw.m_ptr), 0, m_row_size*m_hwm);
    m_hwm = 0ul;
}

void Table_NEW::clear_row(const size_t &irow) {
    std::memset(row_ptr(irow), 0, m_tight_row_size);
}

std::string Table_NEW::to_string(std::string delimiter) const {
    std::string res;
    for (size_t irow=0ul; irow<m_hwm; ++irow) {
        for (const auto field: m_fields) res+=field->to_string(irow)+delimiter+" ";
        res+='\n';
    }
    return res;
}

void Table_NEW::print() const {
    std::cout << to_string("|") << std::endl;
}

std::string Table_NEW::field_details(size_t width) const {
    std::string res;
    for (size_t i = 0ul; i < m_fields.size(); ++i) {
        std::string desc;
        if (!m_fields[i]->m_description.empty()) desc = " (\""+m_fields[i]->m_description+"\")";
        res += "\nField " + std::to_string(i) + desc + ":\n";
        for (auto pair: m_fields[i]->details())
            res += "\t" + utils::padded_string(pair.first, width) + ": " + utils::padded_string(pair.second, width)+"\n";
    }
    return res;
}

void Table_NEW::print_field_details(size_t width) const {
    std::cout << field_details(width) << std::endl;
}
