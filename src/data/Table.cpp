//
// Created by Robert John Anderson on 2020-02-09.
//

#include "Table.h"

#include <iostream>
#include <src/utils.h>

Table::Table(Specification spec, size_t nrow_initial, size_t n_segment,
        float nrow_growth_factor, Table* const shared_table) :
        m_spec(spec.commit()),
        m_row_length(m_spec.m_total_datawords_used),
        m_nsegment(n_segment),
        m_nrow_growth_factor(nrow_growth_factor),
        m_highwatermark(m_nsegment, 0ul),
        m_segment_dataword_offsets(m_nsegment, 0ul),
        m_shared_table(shared_table),
        m_segment_mutex(m_nsegment) {
    for (auto &i: m_segment_mutex) omp_init_lock(&i);
    grow(nrow_initial);
}

Table::~Table() {
    for (auto &i: m_segment_mutex) omp_destroy_lock(&i);
}

void Table::print() const {
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
            for (auto ibitfield{0ul}; ibitfield < m_spec.m_bitfield_lengths.size(); ++ibitfield) {
                std::cout << bitfield_to_string(irow, ibitfield);
            }
            std::cout << std::endl;
        }
    }
}

const std::vector<size_t> &Table::highwatermark() const {
    return m_highwatermark;
}

const std::vector<size_t> &Table::segment_dataword_offsets() const {
    return m_segment_dataword_offsets;
}

size_t Table::row_length() const {
    return m_row_length;
}

size_t Table::get_irow(size_t isegment, size_t isegmentrow, size_t nrow) const {
    return isegment * (nrow ? nrow : m_nrow) + isegmentrow;
}

size_t Table::get_idataword_begin_row(size_t isegment, size_t isegmentrow, size_t nrow) const {
    return get_irow(isegment, isegmentrow, nrow) * m_row_length;
}

size_t Table::get_ndataword_in_segment(size_t isegment) const {
    return m_highwatermark[isegment] * m_row_length;
}

void Table::move_segment(size_t isegment, size_t nrow_old, size_t nrow_new) {
    auto n = get_ndataword_in_segment(isegment);
    auto begin = m_data.begin() + get_idataword_begin_row(isegment, 0, nrow_old);
    std::vector<defs::data_t> tmp(n);
    tmp.assign(begin, begin + n);
    begin = m_data.begin() + get_idataword_begin_row(isegment, 0, nrow_new);
    m_data.insert(begin, tmp.begin(), tmp.end());
}

void Table::grow(size_t nrow_initial) {
    size_t nrow;
    if (nrow_initial) nrow = nrow_initial;
    else nrow = (size_t) (m_nrow * m_nrow_growth_factor);
    m_data.resize(nrow * m_nsegment * m_row_length, 0);
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
        m_segment_dataword_offsets[i] = i * m_row_length * m_nrow;
    }
}

BitfieldNew Table::bitfield_view(const size_t &isegment, const size_t &isegmentrow, const size_t &ientry) const{
    return BitfieldNew(m_spec.m_bitfield_lengths[ientry],
                       (defs::data_t *) dataptr(isegment, isegmentrow, ientry));
}

BitfieldNew Table::bitfield_view(const size_t &irow, const size_t &ientry) const{
    return bitfield_view(0, irow, ientry);
}

size_t Table::grab_rows(size_t isegment, size_t nrow){
    /*
     * advance the highwatermark of the segment without mutex
     */
    size_t first_row;
    first_row = m_highwatermark[isegment];
    auto new_hwm = m_highwatermark[isegment] + nrow;
    if (new_hwm<m_nrow) m_highwatermark[isegment] = new_hwm;
    else {

    }
}

size_t Table::claim_rows(size_t isegment, size_t nrow) {
    size_t first_row;
    lock_acquire(isegment);
    first_row = m_highwatermark[isegment];
    auto new_hwm = m_highwatermark[isegment] + nrow;
    if (new_hwm<m_nrow) m_highwatermark[isegment] = new_hwm;
    else {}//handle_segment_overflow(isegment, first_row);
    lock_release(isegment);
    return first_row;
}

void *Table::row_dataptr(const size_t &isegment, const size_t &isegmentrow) const {
    return (char *) (m_data.data() + get_idataword_begin_row(isegment, isegmentrow));
}

void *Table::row_dataptr(const size_t &irow) const {
    return row_dataptr(0, irow);
}

const Specification &Table::spec() const {
    return m_spec;
}

const size_t Table::nsegment() const {
    return m_nsegment;
}

bool Table::send_to(Table &recv) const {
    MPIWrapper mpi;
    assert(m_nsegment==mpi.nrank());
    assert(recv.nsegment()==1);
    assert(m_spec==recv.spec());

    defs::inds sendcounts = highwatermark();
    for (auto &i : sendcounts) i*=row_length();
    defs::inds recvcounts(mpi.nrank(), 0ul);
    mpi.all_to_all(sendcounts, recvcounts);
    auto &senddispls = segment_dataword_offsets();
    /*
     * place the received data contiguously in the recv buffer:
     */
    defs::inds recvdispls(mpi.nrank(), 0ul);
    for (auto i{1ul}; i<mpi.nrank(); ++i)
        recvdispls[i] = recvdispls[i-1]+sendcounts[i-1];

    return mpi.all_to_allv(baseptr(), sendcounts, senddispls,
                    recv.baseptr(), recvcounts, recvdispls);
}
