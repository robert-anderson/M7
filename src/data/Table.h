//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include "src/defs.h"
#include <omp.h>
#include <cstring>
#include <assert.h>
#include "src/parallel/MPIWrapper.h"
#include "Specification.h"
#include "BitfieldNew.h"
#include "NumericView.h"
#include "src/fermion/Determinant.h"

class Table {
    const Specification m_spec;
    const size_t &m_row_length;
    const size_t m_nsegment;
public:
    const Specification &spec() const;

    const size_t nsegment() const;

private:
    float m_nrow_growth_factor;
    size_t m_nrow{0ul};
    std::vector<defs::data_t> m_data{};
    defs::inds m_highwatermark;
    defs::inds m_segment_dataword_offsets;
    /*
     * if this Table is a thread's private staging space for later inclusion within
     * a shared Table, provide a pointer to that table so the transfer of data rows
     * can be handled when a segment of this table is about to overflow
     */
    Table *const m_shared_table;
    /*
     * When handling overflow events involving private staging Tables, we must ensure
     * that the rest of the team of threads are locked out from:
     *  1. a segment of the table when one thread is claiming transfer space
     *  2. the whole shared table if it is subject to its own overflow, and is
     *     undergoing expansion via the grow method
     *  In the former case, locking the appropriate element of m_segment_mutex is sufficient.
     *  In the latter, a lock must be acquired on every element of that array to prevent all
     *  shared table-related operations until the reallocation is complete
     */
    std::vector<omp_lock_t> m_segment_mutex;

public:
    Table(Specification spec, size_t nrow_initial, size_t n_segment = 1,
          float nrow_growth_factor = 2.0, Table *const shared_table = nullptr);

    ~Table();

    size_t get_irow(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_idataword_begin_row(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_ndataword_in_segment(size_t isegment) const;

    void move_segment(size_t isegment, size_t nrow_old, size_t nrow_new);

    void grow(size_t nrow_initial = 0);

    template<typename T>
    using view_t = typename std::conditional<
            std::is_same<T, BitfieldNew>::value || std::is_same<T, Determinant>::value,
            T, NumericView<T> >::type;

    template<typename T>
    view_t<T> view(const size_t &isegment, const size_t &irow, const size_t &ientry) const {
        if (ientry >= m_spec.m_numeric_lengths[itype<T>])
            return NumericView<T>(nullptr, 0);
        return NumericView<T>(
                (char *) (m_data.data() + (isegment * m_nrow + irow) * m_row_length +
                          m_spec.m_numeric_offsets[itype<T>]) + ientry * sizeof(T),
                m_spec.m_numeric_lengths[itype<T>] - ientry);
    }

    template<typename T>
    view_t<T> view(const size_t &isegment, const size_t &irow) const {
        return view<T>(isegment, irow, 0);
    }

    template<typename T>
    view_t<T> view(const size_t &irow) const {
        return view<T>(0, irow);
    }

    size_t grab_rows(size_t isegment, size_t nrow = 1);

    size_t claim_rows(size_t isegment, size_t nrow = 1);

    size_t request_rows(size_t isegment, size_t nrow = 1);

    size_t row_length() const;

private:

    void lock_acquire(const size_t isegment) {
        omp_set_lock(&m_segment_mutex[isegment]);
    }

    void lock_acquire() {
        /*
         * lock the entire table. Note that this must be undertaken critically
         */
#pragma omp critical
        for (auto isegment{0ul}; isegment < m_nsegment; ++isegment) lock_acquire(isegment);
    }

    void lock_release(const size_t isegment) {
        omp_unset_lock(&m_segment_mutex[isegment]);
    }

    size_t handle_overflow(const size_t &isegment, size_t &row_claimed) {
        /*
         * either this is a shared table
         */
        //if
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

    const defs::data_t *baseptr() const {
        return m_data.data();
    }

    const std::vector<size_t> &highwatermark() const;

    const size_t total_datawords_used() const;

    const std::vector<size_t> &segment_dataword_offsets() const;

    bool send_to(Table &recv) const;

public:
    void print() const;

};

template<>
BitfieldNew Table::view<BitfieldNew>(const size_t &isegment, const size_t &irow, const size_t &ientry) const;

template<>
BitfieldNew Table::view<BitfieldNew>(const size_t &irow, const size_t &ientry) const;

template<>
BitfieldNew Table::view<BitfieldNew>(const size_t &) const = delete;

template<>
Determinant Table::view<Determinant>(const size_t &isegment, const size_t &irow, const size_t &ientry) const;

template<>
Determinant Table::view<Determinant>(const size_t &irow, const size_t &ientry) const;

template<>
Determinant Table::view<Determinant>(const size_t &) const = delete;

#endif //M7_TABLE_H
