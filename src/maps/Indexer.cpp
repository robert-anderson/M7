//
// Created by Robert John Anderson on 2020-01-10.
//

#include "Indexer.h"

defs::inds shape_to_strides(const defs::inds &shape){
    const auto n = shape.size();
    defs::inds strides(n, 0ul);
    strides[n-1] = 1ul;
    for (auto i=2ul; i<=n; i++){
        strides[n-i] = strides[n-i+1]*shape[n-i];
    }
    return strides;
}

size_t Indexer::inds_to_flat(const defs::inds &inds) const {
    /*
     * if the size of the defs::inds is less than ndim,
     * then we assume that the rightmost indices are all zero.
     */
    assert(inds.size()<=ndim_);
    size_t flat = 0ul;
    for (auto i=0ul; i<inds.size(); i++) flat+=inds[i]*m_strides[i];
    return flat;
}

defs::inds Indexer::flat_to_inds(const size_t &flat) const {
    auto work{flat};
    defs::inds inds(ndim_);
    for (auto i=0ul; i<ndim_; i++){
        inds[i] = work/m_strides[i];
        work%=m_strides[i];
    }
    return inds;
}