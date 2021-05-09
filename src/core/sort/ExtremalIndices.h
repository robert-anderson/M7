//
// Created by rja on 29/11/2020.
//

#ifndef M7_EXTREMALINDICES_H
#define M7_EXTREMALINDICES_H


#include "src/defs.h"
#include <functional>
#include "src/core/table/Table.h"

class ExtremalIndices {
    /*
     * inds will be used as a heap, with popped elements being accumulated in
     * order from the back of the vector
     */
    size_t m_hwm = ~0ul;
    defs::inds m_inds;
    typedef std::function<bool(const size_t &, const size_t &)> cmp_t;

    cmp_t m_cmp_fn;
    size_t m_nfound;

public:

    ExtremalIndices(cmp_t cmp_fn);

    const size_t &nfound() const;

    const size_t* begin() const;

    const size_t &operator[](const size_t &ifound) const;

    void reset(size_t hwm);

    void reset(const TableBase &table);

    void find(size_t nfind);

};

#endif //M7_EXTREMALINDICES_H
