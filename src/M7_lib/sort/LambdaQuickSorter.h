//
// Created by Robert J. Anderson on 29/11/2020.
//

#ifndef M7_LAMBDAQUICKSORTER_H
#define M7_LAMBDAQUICKSORTER_H

#include <functional>

#include <M7_lib/defs.h>
#include <M7_lib/table/Table.h>

struct LambdaQuickSorter {

    uintv_t m_inds;
    comparators::index_cmp_fn_t m_cmp_fn;

    LambdaQuickSorter(comparators::index_cmp_fn_t cmp_fn);

    const uint_t& operator[](const uint_t& i) const;

    void preserve_sort(const uint_t &hwm);

    void preserve_sort(const TableBase &table);

    bool is_preserve_sorted(const uint_t &hwm);

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
    void swap(uint_t ii1, uint_t ii2);

    uint_t partition(uint_t iilo, uint_t iihi);

    void qs(uint_t iilo, uint_t iihi);

    void swap(uint_t ii1, uint_t ii2, TableBase &table);

    uint_t partition(uint_t iilo, uint_t iihi, TableBase &table);

    void qs(uint_t iilo, uint_t iihi, TableBase &table);

};


#endif //M7_QUICKSORTER_H
