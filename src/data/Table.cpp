//
// Created by Robert John Anderson on 2020-02-09.
//

#include "Table.h"

#include <iostream>
#include <src/utils.h>


Table::Table(Specification spec, size_t nrow) :
    m_spec(spec) {
    grow(nrow);
}

Table::Table(Specification spec, size_t nrow, defs::data_t *data_external) :
    m_spec(spec) {
    grow(data_external, nrow);
}

template<>
BitfieldNew Table::view<BitfieldNew>(const size_t &irow, const size_t &ientry) const {
    return BitfieldNew(
        m_spec.m_bitfield_lengths[ientry],
        (defs::data_t *) (m_data + irow * row_length() +
                          m_spec.m_total_numeric_datawords_used +
                          m_spec.m_bitfield_offsets[ientry]));
}

template<>
Determinant Table::view<Determinant>(const size_t &irow, const size_t &ientry) const {
    return Determinant(view<BitfieldNew>(irow, ientry), view<BitfieldNew>(irow, ientry + 1));
}

void Table::zero(size_t irow) {
    if (irow != ~0ul) memset(m_data + irow * row_length(), 0, row_length() * sizeof(defs::data_t));
    else memset(m_data, 0, row_length() * m_nrow * sizeof(defs::data_t));
}


void Table::print(size_t nrow) const {
    if (nrow == ~0ul) nrow = m_nrow;
    const size_t padding = 4ul;
    std::cout << "\nnumber of rows: " << m_nrow << std::endl;
    for (size_t irow=0ul; irow < nrow; ++irow) {
        std::cout << view<std::complex<float>>(irow).to_string(padding);
        std::cout << view<std::complex<double>>(irow).to_string(padding);
        std::cout << view<std::complex<long double>>(irow).to_string(padding);
        std::cout << view<float>(irow).to_string(padding);
        std::cout << view<double>(irow).to_string(padding);
        std::cout << view<long double>(irow).to_string(padding);
        std::cout << view<char>(irow).to_string(padding);
        std::cout << view<short int>(irow).to_string(padding);
        std::cout << view<int>(irow).to_string(padding);
        std::cout << view<long int>(irow).to_string(padding);
        std::cout << view<long long int>(irow).to_string(padding);
        std::cout << view<unsigned char>(irow).to_string(padding);
        std::cout << view<unsigned short int>(irow).to_string(padding);
        std::cout << view<unsigned int>(irow).to_string(padding);
        std::cout << view<unsigned long int>(irow).to_string(padding);
        std::cout << view<unsigned long long int>(irow).to_string(padding);
        std::cout << view<bool>(irow).to_string(padding);

        for (size_t ibitfield=0ul; ibitfield < m_spec.m_bitfield_lengths.size(); ++ibitfield) {
            std::cout << view<BitfieldNew>(irow, ibitfield).to_string(padding);
        }
        std::cout << std::endl;
    }
}

void Table::grow(const size_t &nrow) {
    if (m_data && m_data != m_data_internal.data())
        throw std::runtime_error("Cannot reallocate Table which does not own its buffer");
    /*
     * m_data points to m_data_internal, which only needs to be resized
     */
    m_data_internal.resize(nrow * row_length(), 0);
    m_nrow = nrow;
    m_data = m_data_internal.data();
}

void Table::grow(defs::data_t *const new_ptr, size_t nrow) {
    if (m_data && m_data==m_data_internal.data())
        throw std::runtime_error("Cannot move data of a Table which manages its own buffer");
    /*
     * m_data is external, so we trust that the necessary reallocation
     * has already taken place, now we just need to move the data to new pointer
     */
    if (m_data)
        memmove(m_data, new_ptr, high_water_mark() * row_length() * sizeof(defs::data_t));
    m_data = new_ptr;
    m_nrow = nrow;
}

const Specification &Table::spec() const {
    return m_spec;
}

const size_t Table::nrow() const {
    return m_nrow;
}

const size_t Table::row_length() const {
    return m_spec.m_total_datawords_used;
}

const size_t Table::high_water_mark() const {
    return m_nrow;
}

/*
Table::Table(Specification spec, size_t nrow_initial, size_t m_nsegment,
             float nrow_growth_factor, size_t nrow_mutex_blocks) :
    m_spec(spec.compile()),
    m_row_length(m_spec.m_total_datawords_used),
    m_nsegment(m_nsegment),
    m_nrow_growth_factor(nrow_growth_factor),
    m_highwatermark(m_nsegment, 0ul),
    m_segment_dataword_offsets(m_nsegment, 0ul),
    m_nrow_mutex_blocks(nrow_mutex_blocks),
    m_segmentsafe(nrow_mutex_blocks <= m_nsegment),
    m_rowsafe(nrow_mutex_blocks > 0),
    m_segment_mutex(m_segmentsafe ? m_nsegment : 0),
    m_row_mutex(0) {
    grow(nrow_initial);
}

Table::Table(Table *shared, size_t nrow_initial, size_t nrow_mutex_blocks) :
    Table(shared->spec(), nrow_initial, shared->nsegment(), shared->nrow_growth_factor(), nrow_mutex_blocks) {
    m_shared = shared;
}

Table::Table(Table &shared, size_t nrow_initial, size_t nrow_mutex_blocks) :
    Table(&shared, nrow_initial, nrow_mutex_blocks) {}


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
    if (m_rowsafe) m_row_mutex.grow(nrow - m_nrow);
    if (!nrow_initial) {
        for (auto i{m_nsegment}; i > 1ul; --i) move_segment(i, m_nrow, nrow);
    }
    m_nrow = nrow;
    for (size_t i=0ul; i < m_nsegment; ++i) {
        m_segment_dataword_offsets[i] = i * m_row_length * m_nrow;
    }
}


void *Table::row_dataptr(const size_t &isegment, const size_t &isegmentrow) const {
    return (char *) (m_data.data() + get_idataword_begin_row(isegment, isegmentrow));
}

const Specification &Table::spec() const {
    return m_spec;
}

const size_t Table::nsegment() const {
    return m_nsegment;
}

bool Table::send_to(Table &recv) {
    assert(m_nsegment == mpi::nrank());
    assert(recv.nsegment() == 1);
    assert(m_spec == recv.spec());

    defs::inds sendcounts = highwatermark();
    for (auto &i : sendcounts) i *= row_length();
    defs::inds recvcounts(mpi::nrank(), 0ul);
    mpi::all_to_all(sendcounts, recvcounts);
    auto &senddispls = segment_dataword_offsets();
    defs::inds recvdispls(mpi::nrank(), 0ul);
    for (size_t i=1ul; i < mpi::nrank(); ++i)
        recvdispls[i] = recvdispls[i - 1] + sendcounts[i - 1];

    auto tmp = mpi::all_to_allv(baseptr(), sendcounts, senddispls,
                            recv.baseptr(), recvcounts, recvdispls);
    recv.set_highwatermark(0, (recvdispls.back()+recvcounts.back())/row_length());
    zero();
    return tmp;
}

const defs::inds &Table::highwatermark() const {
    return m_highwatermark;
}

void Table::set_highwatermark(const size_t &isegment, const size_t &value) {
    m_highwatermark[isegment] = value;
}

const size_t Table::nrow() const {
    return m_nrow;
}

const size_t Table::nrow_growth_factor() const {
    return m_nrow_growth_factor;
}

void Table::transfer(const size_t &isegment, const size_t n) {
    assert(m_shared != nullptr);
    size_t irow = 0;//m_shared->safe_push(isegment, n);
    memcpy(m_shared->row_dataptr(isegment, irow),
           row_dataptr(isegment, 0),
           n * row_length() * sizeof(defs::data_t));
    m_highwatermark[isegment] = 0;
    // clearing out is not necessary, remove the following after testing
    zero(isegment);
}

size_t Table::push(const size_t &isegment, const size_t &nrow) {
    size_t tmp = m_highwatermark[isegment];
    if (tmp + nrow > m_nrow) {
        if (m_shared != nullptr) transfer(isegment, tmp);
        else throw std::runtime_error("Out of memory to push into Table segment.");
        tmp = 0;
    }
    m_highwatermark[isegment] += nrow;
    return tmp;
}

size_t Table::push(Mutex &mutex, const size_t &nrow) {
    return push(mutex.index(), nrow);
}

const defs::data_t *Table::baseptr() const {
    return m_data.data();
}

defs::data_t *Table::baseptr() {
    return m_data.data();
}

Table::~Table() {
    if (m_shared) {
        for (size_t isegment = 0ul; isegment < m_nsegment; ++isegment) {
            transfer(isegment, m_highwatermark[isegment]);
        }
    }
}
*/
