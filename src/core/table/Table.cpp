//
// Created by Robert John Anderson on 2020-03-26.
//

#include "Table.h"

Table::Table(size_t nsegment) : m_nsegment(nsegment) {}

char *Table::field_begin(const Field *field, const size_t &irow, const size_t isegment) {
    return ((char*)m_data.data())+irow*m_padded_row_size+isegment*m_segment_size+field->m_offset;
}

void Table::expand(size_t delta_rows) {

    /*
     * add more rows to each segment
     */
    m_data.resize(m_data.size() + delta_rows * m_nsegment * m_padded_row_dsize , 0);
    /*
     * move segments backwards to help avoid overlap. std::move will handle overlap correctly if it occurs
     */
    for (size_t isegment = m_nsegment - 1; isegment > 0; --isegment) {
        std::move(
            m_data.begin() + isegment * m_segment_dsize,
            m_data.begin() + (isegment + 1) * m_segment_dsize,
            m_data.begin() + isegment * m_padded_row_dsize*(m_nrow_per_segment + delta_rows)
        );
    }
    increment_nrow_per_segment(delta_rows);
}

size_t Table::irow(const size_t &irow, const size_t &isegment) const {
    return irow + isegment * m_nrow_per_segment;
}

void Table::zero() {
    std::memset(m_data.data(), 0, m_nsegment*m_padded_row_size* sizeof(defs::data_t));
}

size_t Table::add_field(Field *field) {
    if (m_data.size())
        throw std::runtime_error("Cannot add fields: table already in use");
    if (!m_fields.empty() && m_fields.back()->m_type_index!=field->m_type_index){
        // different to last added type, so advance to next
        roundup_row_size();
    }
    size_t offset = m_row_size;
    increment_row_size(field->m_element_size * field->m_nelement);
    m_fields.push_back(field);
    return offset;
}

void Table::update_last_field() {
    auto last_field = m_fields.back();
    m_fields.pop_back();
    update_row_size(last_field->m_offset);
    add_field(last_field);
}

void Table::print() {
    utils::print(m_data);
}

void Table::update_row_size(size_t size) {
    m_row_size = size;
    m_padded_row_size = defs::cache_line_size*integer_utils::divceil(size, defs::cache_line_size);
    assert(m_padded_row_size%sizeof(defs::data_t)==0);
    m_padded_row_dsize = m_padded_row_size/sizeof(defs::data_t);
}

void Table::increment_row_size(size_t delta) {
    update_row_size(m_row_size+delta);
}

void Table::roundup_row_size() {
    update_row_size(integer_utils::divceil(m_row_size, sizeof(defs::data_t))*sizeof(defs::data_t));
}

void Table::update_nrow_per_segment(size_t nrow) {
    m_nrow_per_segment = nrow;
    m_segment_dsize = m_nrow_per_segment*m_padded_row_dsize;
    m_segment_size = m_segment_dsize*sizeof(defs::data_t);
}

void Table::increment_nrow_per_segment(size_t delta) {
    update_nrow_per_segment(m_nrow_per_segment+delta);
}

