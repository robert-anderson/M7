//
// Created by RJA on 18/09/2020.
//

#ifndef M7_NDSPECIFIER_H
#define M7_NDSPECIFIER_H

#include "Indexer.h"

template<typename selector_t, size_t nind>
struct NdSpecifier : NdSpecifierBase<selector_t, nind> {
    // TODO: perfect forwarding
    /*
     * the first indexer has 1 specified index,
     * i.e. nind-1 unspecified indices.
     */
    typedef Indexer<selector_t, nind, nind - 1> first_indexer_t;
    template<typename ...Args>
    NdSpecifier(selector_t selector, const size_t &first, Args ...shape) :
    NdSpecifierBase<selector_t, nind>(selector, first, shape...){}

    first_indexer_t operator[](size_t i) {
        return first_indexer_t(*this, i);
    }
};

template<typename selector_t>
struct NdSpecifier<selector_t, 1ul> : NdSpecifierBase<selector_t, 1ul> {
    // TODO: perfect forwarding
    NdSpecifier(selector_t selector, const size_t &first) :
            NdSpecifierBase<selector_t, 1ul>(selector, first){}

    typedef NdSpecifierBase<selector_t, 1ul> base_t;
    typename selector_t::accessor_t operator[](const size_t& i) {
        return base_t::select(i);
    }
};

// scalar case
template<typename selector_t>
struct NdSpecifier<selector_t, 0ul> : selector_t::accessor_t {
    NdSpecifier<selector_t, 0ul>(): selector_t::accessor_t
    base_t::select(i)
    NdSpecifier<selector_t, 0ul> & operator= ( class_name && )
};




#endif //M7_NDSPECIFIER_H
