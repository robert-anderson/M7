//
// Created by Robert John Anderson on 2020-02-04.
//

#include "DataTable.h"

DataTable::DataTable(
        const std::array<size_t, nnumeric> &numeric_lengths,
        const std::vector<size_t> &bitfield_lengths,
        size_t nrow_initial, size_t n_segment, float nrow_growth_factor) :
        m_numeric_lengths(numeric_lengths),
        m_numeric_offsets(numeric_offsets(numeric_lengths)),
        m_numeric_datawords_used(numeric_datawords_used(numeric_lengths)),
        m_bitfield_lengths(bitfield_lengths),
        m_bitfield_offsets(bitfield_offsets(bitfield_lengths)),
        m_bitfield_datawords_used(bitfield_datawords_used(bitfield_lengths)),
        m_total_numeric_datawords_used(m_numeric_offsets.back() + m_numeric_datawords_used.back()),
        m_total_bitfield_datawords_used(
                bitfield_lengths.empty() ? 0 : m_bitfield_offsets.back() + m_bitfield_datawords_used.back()),
        m_total_datawords_used(m_total_numeric_datawords_used + m_total_bitfield_datawords_used),
        m_nsegment(n_segment),
        m_nrow_growth_factor(nrow_growth_factor),
        m_highwatermark(m_nsegment, 0ul),
        m_segment_dataword_offsets(m_nsegment, 0ul),
        m_segment_mutex(m_nsegment) {
    for (auto &i: m_segment_mutex) omp_init_lock(&i);
    grow(nrow_initial);
}


DataTable::DataTable(const std::array<size_t, nnumeric> &numeric_lengths,
                     size_t nrow_initial, size_t nsegment, float nrow_growth_factor) :
        DataTable(numeric_lengths, defs::inds(), nrow_initial, nsegment, nrow_growth_factor) {}

DataTable::DataTable(const std::vector<size_t> &bitfield_lengths,
                     size_t nrow_initial, size_t nsegment, float nrow_growth_factor) :
        DataTable(std::array<size_t, nnumeric>{}, bitfield_lengths, nrow_initial, nsegment, nrow_growth_factor) {}


DataTable::~DataTable() {
    for (auto &i: m_segment_mutex) omp_destroy_lock(&i);
}

void DataTable::print() const {
    std::cout << "\nnumber of segments: " << m_nsegment;
    std::cout << "\nnumber of rows per segment: " << m_nrow << std::endl;
    for (auto isegment{0ul}; isegment < m_nsegment; ++isegment) {
        std::cout << "\nsegment " << isegment << std::endl;
        for (auto isegmentrow{0ul}; isegmentrow < m_nrow; ++isegmentrow) {
            auto irow = get_irow(isegment, isegmentrow);
            std::cout << numeric_group_to_string<std::complex<float>>(irow);
            std::cout << numeric_group_to_string<std::complex<double>>(irow);
            std::cout << numeric_group_to_string<std::complex<long double>>(irow);
            std::cout << numeric_group_to_string<float>(irow);
            std::cout << numeric_group_to_string<double>(irow);
            std::cout << numeric_group_to_string<long double>(irow);
            std::cout << numeric_group_to_string<char>(irow);
            std::cout << numeric_group_to_string<short int>(irow);
            std::cout << numeric_group_to_string<int>(irow);
            std::cout << numeric_group_to_string<long int>(irow);
            std::cout << numeric_group_to_string<unsigned char>(irow);
            std::cout << numeric_group_to_string<unsigned short int>(irow);
            std::cout << numeric_group_to_string<unsigned int>(irow);
            std::cout << numeric_group_to_string<unsigned long int>(irow);
            for (auto ibitfield{0ul}; ibitfield < m_bitfield_lengths.size(); ++ibitfield) {
                std::cout << bitfield_to_string(irow, ibitfield);
            }
            std::cout << std::endl;
        }
    }
}

const std::vector<size_t> &DataTable::highwatermark() const {
    return m_highwatermark;
}

const std::vector<size_t> &DataTable::segment_dataword_offsets() const {
    return m_segment_dataword_offsets;
}

const size_t DataTable::total_datawords_used() const {
    return m_total_datawords_used;
}

size_t DataTable::get_irow(size_t isegment, size_t isegmentrow, size_t nrow) const {
    return isegment * (nrow ? nrow : m_nrow) + isegmentrow;
}

size_t DataTable::get_idataword_begin_row(size_t isegment, size_t isegmentrow, size_t nrow) const {
    return get_irow(isegment, isegmentrow, nrow) * m_total_datawords_used;
}

size_t DataTable::get_ndataword_in_segment(size_t isegment) const {
    return m_highwatermark[isegment] * m_total_datawords_used;
}

void DataTable::move_segment(size_t isegment, size_t nrow_old, size_t nrow_new) {
    auto n = get_ndataword_in_segment(isegment);
    auto begin = m_data.begin() + get_idataword_begin_row(isegment, 0, nrow_old);
    std::vector<defs::data_t> tmp(n);
    tmp.assign(begin, begin + n);
    begin = m_data.begin() + get_idataword_begin_row(isegment, 0, nrow_new);
    m_data.insert(begin, tmp.begin(), tmp.end());
}

void DataTable::grow(size_t nrow_initial) {
    size_t nrow;
    if (nrow_initial) nrow = nrow_initial;
    else nrow = (size_t) (m_nrow * m_nrow_growth_factor);
    m_data.resize(nrow * m_nsegment * m_total_datawords_used, 0);
    /*
    * move all segments if they number more than one, starting with the one
    * nearest the end of the buffer
    */
    if (!nrow_initial) {
        for (auto i{m_nsegment}; i > 1ul; --i) move_segment(i, m_nrow, nrow);
    }
    m_nrow = nrow;
    /*
     * update the offsets for the beginnings of the segments
     */
    for (auto i{0ul}; i < m_nsegment; ++i) {
        m_segment_dataword_offsets[i] = i * m_total_datawords_used * m_nrow;
    }
}

size_t DataTable::claim_rows(size_t isegment, size_t nrow) {
    size_t first_row;
    omp_set_lock(&m_segment_mutex[isegment]);
    first_row = m_highwatermark[isegment];
    m_highwatermark[isegment] += nrow;
    omp_unset_lock(&m_segment_mutex[isegment]);
    return first_row;
}

void *DataTable::row_dataptr(const size_t &isegment, const size_t &isegmentrow) const {
    return (char *) (m_data.data() + get_idataword_begin_row(isegment, isegmentrow));
}

void *DataTable::row_dataptr(const size_t &irow) const {
    return row_dataptr(0, irow);
}
