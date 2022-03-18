//
// Created by rja on 29/11/2020.
//

#ifndef M7_LAMBDAQUICKSORTER_H
#define M7_LAMBDAQUICKSORTER_H

#include <defs.h>
#include <functional>
#include "table/Table.h"

struct LambdaQuickSorter {

    defs::inds m_inds;
    comparators::index_cmp_fn_t m_cmp_fn;

    LambdaQuickSorter(comparators::index_cmp_fn_t cmp_fn);

    const size_t& operator[](const size_t& i) const;

    void preserve_sort(const size_t &hwm);

    void preserve_sort(const TableBase &table);

    bool is_preserve_sorted(const size_t &hwm);

    bool is_preserve_sorted(const TableBase &table);

    void reorder_sort(TableBase &table);

    /**
     * run through rows and return false if an example is found where two consecutive rows do not conform to the
     * ordering prescribed by m_cmp_fn.
     * e.g. if the sorting criteria is "ascending", and the consecutive values are (1, 2), then this is test passes
     * since the possible criteria which express this ordering (v_i > v_i-1 and v_i >= v_i-1) are both met
     *
     * however, sorts with equal consecutive values must be considered correctly sorted
     *
     *                     CASES
     *  TEST        a < b   a==b    a > b
     *  a < b         1       0       0
     *
     *
     * @param table
     *  table to check for correct ordering
     * @return
     *  true if the rows in the given table are sorted according to m_cmp_fn
     */
    bool is_reorder_sorted(const TableBase &table);

private:
    void swap(size_t ii1, size_t ii2);

    size_t partition(size_t iilo, size_t iihi);

    void qs(size_t iilo, size_t iihi);

    void swap(size_t ii1, size_t ii2, TableBase &table);

    size_t partition(size_t iilo, size_t iihi, TableBase &table);

    void qs(size_t iilo, size_t iihi, TableBase &table);

};


#endif //M7_QUICKSORTER_H
