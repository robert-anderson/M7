//
// Created by rja on 02/10/2020.
//

#include <src/core/util/utils.h>
#include "Table.h"
#include "FieldBase.h"

Table::Table(BufferWindow bw) : m_bw(std::move(bw)){}

Table::Table() : Table(BufferWindow()){}

void Table::move(BufferWindow new_bw) {
    if (m_bw) std::memmove(m_bw.m_ptr, new_bw.m_ptr, sizeof(defs::data_t) * std::min(m_bw.m_dsize, new_bw.m_dsize));
    m_bw = new_bw;
    m_nrow = (sizeof(defs::data_t)*m_bw.m_dsize)/m_row_size;
}

size_t Table::push_back() {
    ASSERT(m_hwm<m_nrow);
    return m_hwm++;
}

void Table::add_field(FieldBase *field) {
    m_ncacheline = integer_utils::divceil(m_tight_row_size, defs::ncacheline_byte);
    m_row_size = m_ncacheline*defs::ncacheline_byte;
    m_fields.push_back(field);
}

char *Table::row_ptr(const size_t &irow) const {
    ASSERT(irow<m_hwm);
    ASSERT(m_bw.m_dsize);
    return (char*)(m_bw.m_ptr)+irow*m_row_size;
}

void Table::clear() {
    std::memset((char*)(m_bw.m_ptr), 0, m_row_size*m_hwm);
    m_hwm = 0ul;
}

void Table::clear_row(const size_t &irow) {
    std::memset(row_ptr(irow), 0, m_tight_row_size);
}

std::string Table::to_string(std::string delimiter) const {
    std::string res;
    for (size_t irow=0ul; irow<m_hwm; ++irow) {
        for (const auto field: m_fields) res+=field->to_string(irow)+delimiter+" ";
        res+='\n';
    }
    return res;
}

void Table::print() const {
    std::cout << to_string("|") << std::endl;
}
