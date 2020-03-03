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
    memset(m_data + irow * row_length(), 0, row_length() * sizeof(defs::data_t));
}

void Table::zero(){
    memset(m_data, 0, row_length() * m_nrow * sizeof(defs::data_t));
}


void Table::print(size_t nrow) const {
    const size_t padding = 4ul;
    if (nrow==m_nrow)
        std::cout << "\nnumber of total rows: " << m_nrow << std::endl;
    else
        std::cout << "\nnumber of rows shown: " << nrow << std::endl;
    for (size_t irow=0ul; irow < nrow; ++irow) {
        std::cout << utils::num_to_string(irow, 6) << " | ";
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

void Table::print() const {
    print(m_nrow);
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

size_t Table::nrow() const {
    return m_nrow;
}

size_t Table::row_length() const {
    return m_spec.m_total_datawords_used;
}

size_t Table::high_water_mark() const {
    return m_nrow;
}
