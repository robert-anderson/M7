//
// Created by Robert J. Anderson on 09/05/2021.
//

#ifndef M7_LOCALEXTREMALROWS_H
#define M7_LOCALEXTREMALROWS_H


#include <M7_lib/field/Fields.h>

#include <utility>
#include "ExtremalIndices.h"

/**
 * Applies the ExtremalIndices class to Table-derived types with numeric fields
 * @tparam T
 *  numerical type of data
 * @tparam nind
 *  number of elements in the shape of the numerical field (0 for scalar)
 */
template<typename T, uint_t nind=0ul>
struct LocalExtremalRows {
    /**
     * references to two order-determining fields in different working rows
     */
    field::Numbers<T, nind> &m_field1;
    field::Numbers<T, nind> &m_field2;
    /**
     * comparator used to determine order of values.
     * if the comparator returns true, then field1 is superior to field2 in the ordering
     */
    const comparators::value_cmp_fn_t<T> m_value_cmp_fn;
    /**
     * the flat elemental indices of the number fields to compare (via arithmetic mean)
     */
    const uintv_t m_inds_to_cmp;
    /**
     * implements the underlying heapsort algorithm
     */
    ExtremalIndices m_xinds;

    LocalExtremalRows(field::Numbers<T, nind> &field1, field::Numbers<T, nind> &field2,
                       bool largest, bool absval, uintv_t inds_to_cmp) :
            m_field1(field1), m_field2(field2),
            m_value_cmp_fn(comparators::get_value_cmp_fn<T>(absval, largest)), m_inds_to_cmp(std::move(inds_to_cmp)),
            m_xinds(comparators::make_num_field_cmp_fn(m_field1, m_field2, m_value_cmp_fn, m_inds_to_cmp)) {
        REQUIRE_NE(field1.m_row, field2.m_row, "work fields rows must point to different rows");
        REQUIRE_EQ(field1.m_row->m_table, field2.m_row->m_table, "both work rows must point to the same table");
        reset();
    }

    LocalExtremalRows(field::Numbers<T, nind> &field1, field::Numbers<T, nind> &field2,
                      bool largest, bool absval, uint_t ind_to_cmp):
            LocalExtremalRows(field1, field2, largest, absval, uintv_t{ind_to_cmp}){}

    void reset() {
        m_xinds.reset(*m_field1.m_row->m_table);
    }

    /**
     * @param nfind
     *  number of local extremal rows to find
     */
    void find(uint_t nfind) {
        m_xinds.find(nfind);
    }

    const uint_t& operator[](uint_t i) const {
        return m_xinds[i];
    }

    uint_t nfound() const {
        return m_xinds.nfound();
    }

    uint_t nremain() const {
        return m_xinds.nremain();
    }

    uint_t nrow_nonzero() const {
        DEBUG_ASSERT_EQ(m_xinds.m_nind, m_field1.m_row->m_table->nrecord(), "inconsistent number of indices");
        return m_xinds.m_nind;
    }

    T get_value(uint_t i) const {
        DEBUG_ASSERT_LT(i, nrow_nonzero(), "row index is OOB wrt number of non-zero rows");
        DEBUG_ASSERT_LT(i, nfound(), "row index is OOB wrt number of found extremal rows");
        m_field1.m_row->jump(m_xinds[i]);
        return m_field1.sum_over(m_inds_to_cmp);
    }

    bool cmp_values(const T& v1, const T& v2) const {
        return m_value_cmp_fn(v1, v2);
    }
};


#endif //M7_LOCALEXTREMALROWS_H
