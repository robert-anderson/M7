//
// Created by Robert John Anderson on 2020-03-26.
//

#include "Table.h"

Table::Table(size_t nsegment) : m_nsegment(nsegment), m_segment_doffsets(nsegment, 0ul) {}

char *Table::field_begin(const Field *field, const size_t &irow, const size_t isegment) {
    ASSERT(is_allocated());
    return ((char *) m_data.data()) + irow * m_padded_row_size + isegment * m_segment_size + field->m_offset;
}

char *Table::row_begin(const size_t &irow, const size_t isegment) {
    return ((char *) m_data.data()) + irow * m_padded_row_size + isegment * m_segment_size;
}

void Table::expand(size_t delta_nrow) {
    /*
     * add more rows to each segment
     */
    m_data.resize((m_nrow_per_segment + delta_nrow) * m_nsegment * m_padded_row_dsize, 0);
    /*
     * move segments backwards to help avoid overlap. std::move will handle overlap correctly if it occurs
     */
    for (size_t isegment = m_nsegment - 1; isegment > 0; --isegment) {
        std::move(
                m_data.begin() + isegment * m_segment_dsize,
                m_data.begin() + (isegment + 1) * m_segment_dsize,
                m_data.begin() + isegment * m_padded_row_dsize * (m_nrow_per_segment + delta_nrow)
        );
    }
    increment_nrow_per_segment(delta_nrow);
    for (size_t isegment = 1ul; isegment < m_nsegment; ++isegment) {
        m_segment_doffsets[isegment] = m_segment_doffsets[isegment - 1] + m_padded_row_dsize * m_nrow_per_segment;
    }
}

void Table::resize(size_t nrow) {
    if (nrow > m_nrow_per_segment) expand(nrow - m_nrow_per_segment);
}

size_t Table::irow(const size_t &irow, const size_t &isegment) const {
    return irow + isegment * m_nrow_per_segment;
}

void Table::zero() {
    std::memset(m_data.data(), 0, m_nsegment * m_nrow_per_segment * m_padded_row_size);
}

void Table::zero_row(const size_t &irow, const size_t &isegment) {
    std::memset(row_begin(irow, isegment), 0, m_padded_row_size);
}

size_t Table::add_field(Field *field) {
    if (m_data.size())
        throw std::runtime_error("Cannot add fields: table already in use");
    if (!m_fields.empty() && m_fields.back()->m_type_index != field->m_type_index) {
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
    std::cout << to_string() << std::endl;
}

const size_t &Table::nrow_per_segment() const {
    return m_nrow_per_segment;
}

void Table::update_row_size(size_t size) {
    m_row_size = size;
    m_padded_row_size = defs::cache_line_size * integer_utils::divceil(size, defs::cache_line_size);
    ASSERT(m_padded_row_size % sizeof(defs::data_t) == 0);
    m_padded_row_dsize = m_padded_row_size / sizeof(defs::data_t);
}

void Table::increment_row_size(size_t delta) {
    update_row_size(m_row_size + delta);
}

void Table::roundup_row_size() {
    update_row_size(integer_utils::divceil(m_row_size, sizeof(defs::data_t)) * sizeof(defs::data_t));
}

void Table::update_nrow_per_segment(size_t nrow) {
    m_nrow_per_segment = nrow;
    m_segment_dsize = m_nrow_per_segment * m_padded_row_dsize;
    m_segment_size = m_segment_dsize * sizeof(defs::data_t);
}

void Table::increment_nrow_per_segment(size_t delta) {
    update_nrow_per_segment(m_nrow_per_segment + delta);
}

bool Table::compatible_with(const Table &other) const {
    if (m_fields.size() != other.m_fields.size()) return false;
    for (size_t ifield = 0ul; ifield < m_fields.size(); ++ifield) {
        if (!m_fields[ifield]->compatible_with(*other.m_fields[ifield])) return false;
    }
    return m_row_size == other.m_row_size &&
           m_padded_row_size == other.m_padded_row_size &&
           m_padded_row_dsize == other.m_padded_row_dsize;
}

bool Table::is_allocated() const { return m_data.data() != nullptr; }

std::string Table::row_to_string(size_t irow, size_t isegment) {
    std::string result;
    for (auto &field : m_fields) {
        for (size_t ielement = 0ul; ielement < field->m_nelement; ++ielement) {
            result += field->to_string(irow, isegment, ielement) + "  ";
        }
    }
    return result;
}

std::string Table::to_string(const defs::inds &nrows) {
    std::string result = "\nTABLE\n";
    result += "# fields: " + std::to_string(m_fields.size()) + "\n";
    for (size_t ifield = 0ul; ifield < m_fields.size(); ++ifield) {
        result += "  " + std::to_string(ifield) + ": " + m_fields[ifield]->description() + "\n";
    }
    for (size_t isegment = 0ul; isegment < m_nsegment; ++isegment) {
        if (!m_nsegment) result += "SEGMENT " + std::to_string(isegment) + "\n";
        result += "# rows: " + std::to_string(nrows[isegment]) + "\n";
        for (size_t irow = 0ul; irow < nrows[isegment]; ++irow) {
            result += utils::num_to_string(unsigned(irow)) + ":  " + row_to_string(irow, isegment) + "\n";
        }
    }
    return result;
}

std::string Table::to_string() {
    return to_string(defs::inds(m_nsegment, m_nrow_per_segment));
}

void Table::print_row(size_t irow, size_t isegment) {
    std::cout << row_to_string(irow, isegment) << std::endl;
}

size_t Table::dsize() const {
    return m_segment_dsize * m_nsegment;
}
