//
// Created by rja on 09/05/2021.
//

#ifndef M7_LOCALEXTREMALROWS_H
#define M7_LOCALEXTREMALROWS_H


#include <src/core/util/utils.h>
#include <src/core/field/Fields.h>
#include "ExtremalIndices.h"

/**
 * Applies the ExtremalIndices class to Table-derived types with numeric fields
 * @tparam row_t
 *  row template argument
 * @tparam T
 *  numerical type of data
 * @tparam nind
 *  number of elements in the shape of the numerical field (0 for scalar)
 */
template<typename row_t, typename T, size_t nind=0>
struct LocalExtremalRows {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    /**
     * once-initialized, frequently assigned working row
     */
    row_t m_work_row;
    /**
     * reference to order-determining field in m_work_row;
     */
    fields::Numbers <T, nind> &m_work_row_field;
    /**
     * once-initialized, frequently assigned working row for comparisons
     */
    row_t m_work_row_cmp;
    /**
     * reference to order-determining field in m_work_row_cmp;
     */
    fields::Numbers <T, nind> &m_work_row_cmp_field;
    /**
     * get the largest values in the local table
     */
    const bool m_largest;
    /**
     * take the std::abs of numeric data before comparison
     */
    const bool m_absval;
    /**
     * the element index of the number field to compare.
     * if this is ~0ul, the comparator sorts based on all columns lexically, with 0 being the major element index
     */
    const size_t m_ielement_cmp;
    /**
     * comparator used to determine order of values
     */
    const std::function<bool(const T &, const T &)> m_value_cmp_fn;
    /**
     * implements the underlying heapsort algorithm
     */
    ExtremalIndices m_xinds;

    LocalExtremalRows(row_t &row, fields::Numbers<T, nind> &field, bool largest, bool absval, size_t ielement_cmp=~0ul) :
            m_work_row(row),
            m_work_row_field(fields::identify(m_work_row, row, field)),
            m_work_row_cmp(row),
            m_work_row_cmp_field(fields::identify(m_work_row_cmp, row, field)),
            m_largest(largest), m_absval(absval), m_ielement_cmp(ielement_cmp),
            m_value_cmp_fn(comparators::make_value_cmp_fn<T>(absval, largest)),
            m_xinds(comparators::make_num_field_row_cmp_fn(
                    m_work_row, m_work_row_field,
                    m_work_row_cmp, m_work_row_cmp_field, m_value_cmp_fn, ielement_cmp)) {
        ASSERT(static_cast<const FieldBase&>(m_work_row_field).belongs_to_row(m_work_row));
        ASSERT(static_cast<const FieldBase&>(m_work_row_cmp_field).belongs_to_row(m_work_row_cmp));
        reset();
    }

    void reset() {
        m_xinds.reset(static_cast<const Row &>(m_work_row).m_table->m_hwm);
    }

    /**
     * @param nfind
     *  number of local extremal rows to find
     */
    void find(size_t nfind) {
        m_xinds.find(nfind);
    }

    const size_t& operator[](const size_t& i) const {
        return m_xinds[i];
    }

    const size_t &nfound() const {
        return m_xinds.nfound();
    }

    size_t nremain() const {
        return m_xinds.nremain();
    }

    size_t hwm() const {
        return static_cast<const Row &>(m_work_row).m_table->m_hwm;
    }


};


#endif //M7_LOCALEXTREMALROWS_H
