//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <src/defs.h>
#include <omp.h>
#include <cstring>
#include "Specification.h"
#include "BitfieldNew.h"


class Table {
    const Specification m_spec;
    const size_t &m_row_length;
    const size_t m_nsegment;
    float m_nrow_growth_factor;
    size_t m_nrow{0ul};
    std::vector<defs::data_t> m_data{};
    defs::inds m_highwatermark;
    defs::inds m_segment_dataword_offsets;
    std::vector<omp_lock_t> m_segment_mutex;

public:
    Table(Specification spec, size_t nrow_initial, size_t n_segment = 1, float nrow_growth_factor = 2.0);

    ~Table();

    size_t get_irow(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_idataword_begin_row(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_ndataword_in_segment(size_t isegment) const;

    void move_segment(size_t isegment, size_t nrow_old, size_t nrow_new);

    void grow(size_t nrow_initial = 0);

    template<typename T>
    void encode(const T &v, const size_t &isegment, const size_t &isegmentrow, const size_t &ientry) {
        assert(ientry < m_spec.m_numeric_lengths[itype<T>]);
        memcpy(dataptr<T>(isegment, isegmentrow, ientry), &v, sizeof(T));
    }

    template<typename T>
    void encode(const T &v, const size_t &irow, const size_t &ientry) {
        encode(v, 0, irow, ientry);
    }

    template<typename T>
    void encode(const std::vector<T> &v, const size_t &isegment, const size_t &isegmentrow) {
        auto length = std::min(v.size(), m_spec.m_numeric_lengths[itype<T>]);
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

    BitfieldNew bitfield_view(const size_t &isegment, const size_t &isegmentrow, const size_t &ientry) const;

    BitfieldNew bitfield_view(const size_t &irow, const size_t &ientry) const;

    template<typename T>
    void decode(T &v, const size_t &irow, const size_t &ientry) const {
        auto offset = m_spec.m_numeric_offsets[itype<T>];
        memcpy(&v, dataptr<T>(irow, ientry), sizeof(T));
    }

    template<typename T>
    void decode(std::vector<T> &v, const size_t &isegment, const size_t &isegmentrow) const {
        auto length = std::min(v.size(), m_spec.m_numeric_lengths[itype<T>]);
        memcpy(v.data(), dataptr<T>(isegment, isegmentrow, 0), length * sizeof(T));
    }

    template<typename T>
    void decode(std::vector<T> &v, const size_t &irow) const {
        decode(v, 0, irow);
    }

    size_t claim_rows(size_t isegment, size_t nrow = 1);

    size_t row_length() const;

private:

    template<typename T>
    void *dataptr(const size_t &isegment, const size_t &isegmentrow, const size_t &ientry) const {
        // pointer to specific entry of numeric type
        return (char *) (m_data.data() + (isegment * m_nrow + isegmentrow) * m_row_length +
                         m_spec.m_numeric_offsets[itype<T>]) + ientry * sizeof(T);
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
                         m_spec.m_total_numeric_datawords_used + m_spec.m_bitfield_offsets[ibitfield]);
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
        for (auto i{0ul}; i < m_spec.m_numeric_lengths[itype<T>]; ++i) {
            decode(tmp, irow, i);
            out += entry_to_string(tmp, padding);
        }
        return out;
    }

    std::string bitfield_to_string(size_t irow, size_t ibitfield, size_t padding = 2) const {
        BitfieldNew bitfield = bitfield_view(irow, ibitfield);
        auto tmp_string = bitfield.to_string();
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }

public:
    void print() const;

};

#endif //M7_TABLE_H
