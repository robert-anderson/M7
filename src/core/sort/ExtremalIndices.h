//
// Created by rja on 29/11/2020.
//

#ifndef M7_EXTREMALINDICES_H
#define M7_EXTREMALINDICES_H


#include "src/defs.h"
#include <functional>
#include <src/core/parallel/MPIAssert.h>
#include "src/core/table/TableBase.h"
#include "Comparators.h"

/**
 * given a comparator function between two elements, a total number of elements and a number of elements to find, this
 * class, via the find method will store that number of element indices in the order determined by the sense of the
 * inequality in the comparator
 */
class ExtremalIndices {
    /**
     * the "high water mark" or total number of elements from which to find the indices corresponding to the ordered
     * extremal values
     */
    size_t m_hwm = ~0ul;
    /**
     * after a reset, this is ordered and consecutive. After a find, the first m_nfound are the extremal indices in the
     * requested order
     */
    defs::inds m_inds;
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

public:

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
     * start from scratch with an ordered inds array
     * @param hwm
     *  new high water mark (total number of elements to heap sort a subset from)
     */
    void reset(size_t hwm);
    /**
     * @param table
     *  table object from which the m_hwm is accessed
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
