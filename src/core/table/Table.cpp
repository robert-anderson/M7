//
// Created by rja on 21/10/2020.
//

#include "Table.h"


bool Table::is_full() const {
    return m_hwm==m_nrow;
}

size_t Table::push_back(size_t nrow) {
    if (m_hwm>=m_nrow) expand();
    auto tmp = m_hwm;
    m_hwm+=nrow;
    return tmp;
}

defs::data_t *Table::dbegin() {
    return m_bw.dbegin();
}

const defs::data_t *Table::dbegin() const {
    return m_bw.dbegin();
}

char *Table::begin() {
    return (char*) m_bw.dbegin();
}

const char *Table::begin() const {
    return (const char*)m_bw.dbegin();
}

defs::data_t *Table::dbegin(const size_t &irow) {
    ASSERT(irow<m_hwm)
    ASSERT(m_bw.dbegin())
    return m_bw.dbegin() + irow * m_row_dsize;
}

const defs::data_t *Table::dbegin(const size_t &irow) const {
    ASSERT(irow<m_hwm)
    ASSERT(m_bw.dbegin())
    return m_bw.dbegin() + irow * m_row_dsize;
}

char *Table::begin(const size_t &irow) {
    return (char*)dbegin(irow);
}

const char *Table::begin(const size_t &irow) const {
    return (const char*)dbegin(irow);
}

size_t Table::add_field(const TableField *field) {
    // returns the offset in bytes for the field being added
    auto offset = 0ul;
    if(!m_fields.empty()){
        offset = m_fields.back()->m_offset+m_fields.back()->m_size;
        if (!m_fields.back()->is_same_type_as(*field)){
            // go to next whole dataword
            offset = integer_utils::divceil(offset, defs::nbyte_data)*defs::nbyte_data;
        }
    }

    m_current_byte_offset = offset + field->m_size;
    m_row_dsize = integer_utils::divceil(m_current_byte_offset, defs::nbyte_data);
    m_row_size = m_row_dsize*defs::nbyte_data;

    m_fields.push_back(field);
    return offset;
}

//void Table::move(BufferWindow new_bw) {
//    if (m_bw) std::memmove(new_bw.m_ptr, m_bw.m_ptr, sizeof(defs::data_t) * std::min(m_bw.m_dsize, new_bw.m_dsize));
//    m_bw = new_bw;
//    if (!m_row_dsize) return;
//    m_nrow = m_bw.m_dsize / m_row_dsize;
//}

void Table::clear() {
    std::memset(begin(), 0, m_row_size * m_hwm);
    m_hwm = 0ul;
}

void Table::clear(const size_t &irow) {
    std::memset(begin(irow), 0, m_row_size);
}

bool Table::is_cleared() const {
    return std::all_of(dbegin(), dbegin()+m_row_dsize*m_nrow, [](const defs::data_t& i){return i==0;});
}

bool Table::is_cleared(const size_t &irow) const {
    return std::all_of(dbegin(irow), dbegin(irow)+m_row_dsize, [](const defs::data_t& i){return i==0;});
}

std::string Table::field_details(size_t width) const {
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

void Table::print_field_details(size_t width) const {
    std::cout << field_details(width) << std::endl;
}

size_t Table::bw_dsize() const {
    return m_bw.dsize();
}

void Table::print_contents(const defs::inds *ordering) const {
    const auto n = ordering ? std::min(ordering->size(), m_hwm) : m_hwm;
    for (size_t iirow=0ul; iirow<n; ++iirow){
        auto irow = ordering ? (*ordering)[iirow] : iirow;
        std::cout << irow << ". ";
        for (auto field: m_fields){
            std::cout << field->to_string(irow)+" ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void Table::print_contents(const ExtremalValues &xv) const {
    defs::inds tmp;
    tmp.reserve(xv.nfound());
    for (size_t i=0ul; i<xv.nfound(); ++i) tmp.push_back(xv[i]);
    print_contents(&tmp);
}