//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_DATATABLE_H
#define M7_DATATABLE_H


#include <cstddef>
#include <vector>
#include <array>
#include <complex>
#include <assert.h>
#include <iostream>
#include <cstring>
#include <omp.h>
#include "BitfieldNew.h"

template<typename T>
constexpr size_t type_number = ~0ul;

// complex numbers
template<>
constexpr size_t type_number<std::complex<float>> = 0;
template<>
constexpr size_t type_number<std::complex<double>> = 1;
template<>
constexpr size_t type_number<std::complex<long double>> = 2;


// real numbers
template<>
constexpr size_t type_number<float> = 3;
template<>
constexpr size_t type_number<double> = 4;
template<>
constexpr size_t type_number<long double> = 5;

// signed integers
template<>
constexpr size_t type_number<char> = 6;
template<>
constexpr size_t type_number<short int> = 7;
template<>
constexpr size_t type_number<int> = 8;
template<>
constexpr size_t type_number<long int> = 9;
template<>
constexpr size_t type_number<long long int> = 10;

// unsigned integers
template<>
constexpr size_t type_number<unsigned char> = 11;
template<>
constexpr size_t type_number<unsigned short int> = 12;
template<>
constexpr size_t type_number<unsigned int> = 13;
template<>
constexpr size_t type_number<unsigned long int> = 14;
template<>
constexpr size_t type_number<unsigned long long int> = 15;

// boolean
template<>
constexpr size_t type_number<bool> = 16;


constexpr size_t nnumeric = 17;

constexpr std::array<size_t, nnumeric> type_sizes = {
        sizeof(std::complex<float>),
        sizeof(std::complex<double>),
        sizeof(std::complex<long double>),
        sizeof(float),
        sizeof(double),
        sizeof(long double),
        sizeof(char),
        sizeof(short int),
        sizeof(int),
        sizeof(long int),
        sizeof(long long int),
        sizeof(unsigned char),
        sizeof(unsigned short int),
        sizeof(unsigned int),
        sizeof(unsigned long int),
        sizeof(unsigned long long int),
        sizeof(bool)
};

class DataTable {
    const std::array<size_t, nnumeric> m_numeric_lengths;
    const std::array<size_t, nnumeric> m_numeric_datawords_used;
    const std::array<size_t, nnumeric> m_numeric_offsets;
    const std::vector<size_t> m_bitfield_lengths;
    const std::vector<size_t> m_bitfield_datawords_used;
    const std::vector<size_t> m_bitfield_offsets;
    const size_t m_total_numeric_datawords_used;
    const size_t m_total_bitfield_datawords_used;
    const size_t m_total_datawords_used;

private:
    const size_t m_nsegment;
    float m_nrow_growth_factor;
    size_t m_nrow{0ul};
    std::vector<defs::data_t> m_data{};
    defs::inds m_highwatermark;
    defs::inds m_segment_dataword_offsets;
    std::vector<omp_lock_t> m_segment_mutex;

public:
    DataTable(
            const std::array<size_t, nnumeric> &numeric_lengths,
            const std::vector<size_t> &bitfield_lengths,
            size_t nrow_initial, size_t n_segment = 1, float nrow_growth_factor = 2.0
    );
    DataTable(
            const std::array<size_t, nnumeric> &numeric_lengths,
            size_t nrow_initial, size_t nsegment = 1, float nrow_growth_factor = 2.0
    );
    DataTable(
            const std::vector<size_t> &bitfield_lengths,
            size_t nrow_initial, size_t n_segment = 1, float nrow_growth_factor = 2.0
    );

    ~DataTable();

    size_t get_irow(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_idataword_begin_row(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_ndataword_in_segment(size_t isegment) const;

    void move_segment(size_t isegment, size_t nrow_old, size_t nrow_new);

    void grow(size_t nrow_initial = 0);

    template<typename T>
    void encode(const T &v, const size_t &isegment, const size_t &isegmentrow, const size_t &ientry) {
        assert(ientry < m_numeric_lengths[type_number<T>]);
        memcpy(dataptr<T>(isegment, isegmentrow, ientry), &v, sizeof(T));
    }

    template<typename T>
    void encode(const T &v, const size_t &irow, const size_t &ientry) {
        encode(v, 0, irow, ientry);
    }

    template<typename T>
    void encode(const std::vector<T> &v, const size_t &isegment, const size_t &isegmentrow) {
        auto length = std::min(v.size(), m_numeric_lengths[type_number<T>]);
        memcpy(dataptr<T>(isegment, isegmentrow, 0), v.data(), length * sizeof(T));
    }

    template<typename T>
    void encode(const std::vector<T> &v, const size_t &irow) {
        encode(v, 0, irow);
    }

    template<typename T>
    T *view(const size_t &isegment, const size_t &irow, const size_t &ientry) {
        return (T *) dataptr<T>(isegment, irow, ientry);
    }

    template<typename T>
    T *view(const size_t &irow, const size_t &ientry) {
        return view<T>(0, irow, ientry);
    }

    BitfieldNew view(const size_t &isegment, const size_t &irow, const size_t &ientry) const {
        return BitfieldNew(m_bitfield_lengths[ientry],
                           (defs::data_t *) dataptr(isegment, irow, ientry));
    }

    BitfieldNew view(const size_t &irow, const size_t &ientry) const {
        return view(0, irow, ientry);
    }

    template<typename T>
    void decode(T &v, const size_t &irow, const size_t &ientry) const {
        auto offset = m_numeric_offsets[type_number<T>];
        memcpy(&v, dataptr<T>(irow, ientry), sizeof(T));
    }

    template<typename T>
    void decode(std::vector<T> &v, const size_t &isegment, const size_t &isegmentrow) const {
        auto length = std::min(v.size(), m_numeric_lengths[type_number<T>]);
        memcpy(v.data(), dataptr<T>(isegment, isegmentrow, 0), length * sizeof(T));
    }

    template<typename T>
    void decode(std::vector<T> &v, const size_t &irow) const {
        decode(v, 0, irow);
    }

    size_t claim_rows(size_t isegment, size_t nrow = 1);

private:

    template<typename T>
    void *dataptr(const size_t &isegment, const size_t &isegmentrow, const size_t &ientry) const {
        // pointer to specific entry of numeric type
        return (char *) (m_data.data() + (isegment * m_nrow + isegmentrow) * m_total_datawords_used +
                         m_numeric_offsets[type_number<T>]) + ientry * sizeof(T);
    }

    template<typename T>
    void *dataptr(const size_t &irow, const size_t &ientry) const {
        return dataptr<T>(0, irow, ientry);
    }

    template<typename T>
    void *dataptr(const size_t &irow) const {
        // pointer to beginning of numeric type
        return dataptr<T>(irow, 0);
    }

    void *dataptr(const size_t &isegment, const size_t &isegmentrow, const size_t &ibitfield) const {
        // pointer to specific bitfield
        return (char *) (m_data.data() + get_idataword_begin_row(isegment, isegmentrow) +
                         m_total_numeric_datawords_used + m_bitfield_offsets[ibitfield]);
    }

    void *dataptr(const size_t &irow, const size_t &ibitfield) const {
        return dataptr(0, irow, ibitfield);
    }

    void *row_dataptr(const size_t &isegment, const size_t &isegmentrow) const;

    void *row_dataptr(const size_t &irow) const;

public:
    /*
     * these methods must be exposed for the purpose of MPI communication
     */
    defs::data_t *baseptr() {
        return m_data.data();
    }

    const std::vector<size_t> &highwatermark() const;

    const size_t total_datawords_used() const;

    const std::vector<size_t> &segment_dataword_offsets() const;

private:
    template<typename T>
    std::string entry_to_string(T &entry, size_t padding = 0) const {
        auto tmp_string = std::to_string(entry);
        auto decimal_length = std::numeric_limits<T>::digits10;
        tmp_string.insert(tmp_string.begin(), padding + decimal_length - tmp_string.size(), ' ');
        return tmp_string;
    }

    template<typename T>
    std::string entry_to_string(std::complex<T> &entry, size_t padding = 0) const {
        auto tmp_string = std::to_string(entry.real()) +
                          (entry.imag() < 0 ? "" : "+") + std::to_string(entry.imag()) + "i";
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }

    template<typename T>
    std::string numeric_group_to_string(size_t irow, size_t padding = 4) const {
        std::string out{};
        T tmp;
        for (auto i{0ul}; i < m_numeric_lengths[type_number<T>]; ++i) {
            decode(tmp, irow, i);
            out += entry_to_string(tmp, padding);
        }
        return out;
    }

    std::string bitfield_to_string(size_t irow, size_t ibitfield, size_t padding = 2) const {
        BitfieldNew bitfield = view(irow, ibitfield);
        auto tmp_string = bitfield.to_string();
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }

    static std::array<size_t, nnumeric> numeric_datawords_used(const std::array<size_t, nnumeric> &lengths) {
        std::array<size_t, nnumeric> out{};
        for (auto i{0ul}; i < nnumeric; ++i) {
            out[i] = (lengths[i] * type_sizes[i]) / sizeof(defs::data_t) +
                     (((lengths[i] * type_sizes[i]) % sizeof(defs::data_t)) != 0);
        }
        return out;
    }

    static std::array<size_t, nnumeric> numeric_offsets(const std::array<size_t, nnumeric> &lengths) {
        std::array<size_t, nnumeric> out{};
        auto tmp = numeric_datawords_used(lengths);
        out[0] = 0;
        for (auto i{1ul}; i < nnumeric; ++i) {
            out[i] = out[i - 1] + tmp[i - 1];
        }
        return out;
    }

    static std::vector<size_t> bitfield_datawords_used(const std::vector<size_t> &lengths) {
        std::vector<size_t> out(lengths.size(), 0ul);
        for (auto i{0ul}; i < lengths.size(); ++i) {
            out[i] = lengths[i] / (sizeof(defs::data_t) * 8) + ((lengths[i] % (sizeof(defs::data_t) * 8)) != 0);
        }
        return out;
    }

    static std::vector<size_t> bitfield_offsets(const std::vector<size_t> &lengths) {
        std::vector<size_t> out(lengths.size(), 0ul);
        auto tmp = bitfield_datawords_used(lengths);
        for (auto i{1ul}; i < lengths.size(); ++i) {
            out[i] = out[i - 1] + tmp[i - 1];
        }
        return out;
    }

public:
    void print() const;

};


#endif //M7_DATATABLE_H
