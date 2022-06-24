//
// Created by Robert J. Anderson on 29/11/2020.
//

#ifndef M7_EXTREMALINDICES_H
#define M7_EXTREMALINDICES_H


#include <functional>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/table/TableBase.h>

#include "Comparators.h"

/**
 * given a comparator function between two elements, a total number of elements and a number of elements to find, this
 * class, via the find method will store that number of element indices in the order determined by the sense of the
 * inequality in the comparator
 */
struct ExtremalIndices {
    /**
     * total number of elements from which to find the indices corresponding to the ordered extremal values, for tables
     * this is (the high water mark) - (the number of row on the free stack)
     */
    size_t m_nind = ~0ul;
    /**
     * after a reset, this is ordered and consecutive. After a find, the first m_nfound are the extremal indices in the
     * requested order
     */
    defs::uintv_t m_inds;
    /**
     * the comparator - the lambda technique is powerful because the calling scope can declare a lambda which captures
     * any working objects required and to evaluate the relative order of two elements without this class needing to be
     * concerned with those specifics
     */
    comparators::index_cmp_fn_t m_cmp_fn;
    /**
     * total number of extremal valued elements found so far. never greater than m_hwm.
     */
    size_t m_nfound;

    explicit ExtremalIndices(comparators::index_cmp_fn_t cmp_fn);

    /**
     * @return
     *  the total number of elements found so far over all calls to find
     */
    const size_t &nfound() const;
    /**
     * @return
     *  the number of elements as yet unfound
     */
    size_t nremain() const;
    /**
     * @return
     *  pointer to the last found element index
     */
    const size_t* begin() const;
    /**
     * @param ifound
     *  index of the found element indices in order determined by the comparator
     * @return
     *  the ifound-th found element index
     */
    const size_t &operator[](const size_t &ifound) const;
    /**
     * start from scratch with an ordered uintv_t array
     * @param hwm
     *  new high water mark (total number of elements to sort a subset from)
     * @param inds_ignore
     *  unsorted vector of indices which should be ignored in the partial sort
     */
    void reset(size_t hwm, defs::uintv_t inds_ignore={});
    /**
     * @param table
     *  table object from which the m_hwm and free row indices are accessed
     */
    void reset(const TableBase &table);
    /**
     * execute a partial sort in order to find an additional nfind elements
     * @param nfind
     *  number of additional elements to find
     */
    void find(size_t nfind);
};

#endif //M7_EXTREMALINDICES_H
