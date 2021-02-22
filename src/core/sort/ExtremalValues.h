//
// Created by rja on 29/11/2020.
//

#ifndef M7_EXTREMALVALUES_H
#define M7_EXTREMALVALUES_H


#include "src/defs.h"
#include <functional>
#include "src/core/fieldz/TableZ.h"

class ExtremalValues {
    /*
     * inds will be used as a heap, with popped elements being accumulated in
     * order from the back of the vector
     */
    size_t m_hwm = ~0ul;
    defs::inds m_inds;
    typedef std::function<bool(const size_t &, const size_t &)> comp_t;

    comp_t m_comp_fn;
    size_t m_nfound;

public:

    ExtremalValues(comp_t comp_fn);

    const size_t &nfound() const;

    const size_t* begin() const;

    const size_t &operator[](const size_t &ifound) const;

    void reset(size_t hwm);

    void reset(const TableBaseZ &table);

    void find(size_t nfind);

    const defs::inds& inds() const;

};

#endif //M7_EXTREMALVALUES_H
