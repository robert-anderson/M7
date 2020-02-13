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
#include "MutexVector.h"

class Table {
    const Specification m_spec;
    const size_t &m_row_length;
protected:
    const size_t m_nsegment;
    float m_nrow_growth_factor;
    size_t m_nrow{0ul};
    std::vector<defs::data_t> m_data{};
    defs::inds m_highwatermark;
    defs::inds m_segment_dataword_offsets;

    const size_t m_nrow_mutex_blocks;
    const bool m_segmentsafe;
    const bool m_rowsafe;
    MutexVector m_segment_mutex;
    MutexVector m_row_mutex;

    Table* m_shared = nullptr;

public:

    const Specification &spec() const;

    const size_t nsegment() const;

    const size_t nrow() const;

    const size_t nrow_growth_factor() const;

    Table(Specification spec, size_t nrow_initial, size_t m_nsegment = 1,
          float nrow_growth_factor = 2.0, size_t nrow_mutex_blocks=0);

    Table(Table *shared, size_t nrow_initial, size_t nrow_mutex_blocks=0);
    Table(Table &shared, size_t nrow_initial, size_t nrow_mutex_blocks=0);

    virtual ~Table();

    size_t get_irow(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_idataword_begin_row(size_t isegment, size_t isegmentrow, size_t nrow = 0) const;

    size_t get_ndataword_in_segment(size_t isegment) const;

    void move_segment(size_t isegment, size_t nrow_old, size_t nrow_new);

    void grow(size_t nrow_initial = 0);

    void zero(size_t isegment=~0ul, size_t irow=~0ul);

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

    const defs::inds& highwatermark() const;

    void transfer(const size_t &isegment, const size_t n);

    virtual size_t push(const size_t &isegment, const size_t &nrow);

    virtual size_t safe_push(const size_t &isegment, const size_t &nrow);

    size_t row_length() const;

    void *row_dataptr(const size_t &isegment, const size_t &isegmentrow) const;

    void *row_dataptr(const size_t &irow) const;

    /*
     * these methods must be exposed for the purpose of MPI communication
     */
    defs::data_t *baseptr();

    const defs::data_t *baseptr() const;

    const size_t total_datawords_used() const;

    const std::vector<size_t> &segment_dataword_offsets() const;

    bool send_to(Table &recv) const;

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
