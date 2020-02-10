//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_INDEXER_H
#define M7_INDEXER_H

#include <vector>
#include <iostream>
#include <assert.h>
#include "../defs.h"

defs::inds shape_to_strides(const defs::inds &shape);

class Indexer {
    const defs::inds m_shape, m_strides;
public:
    const size_t ndim_;
    const size_t nelem_;
    Indexer(const defs::inds &shape, size_t ndim):
            m_shape{(assert(ndim<=shape.size()), defs::inds(shape.begin(), shape.end()+(ndim-shape.size())))},
            m_strides{shape_to_strides(m_shape)},
            ndim_{ndim}, nelem_{m_strides[0]*m_shape[ndim-1]}{}
    Indexer(const defs::inds &shape): Indexer(shape, shape.size()){};
    Indexer(const size_t &extent, const size_t &ndim): Indexer(defs::inds(ndim, extent)){};
    size_t inds_to_flat(const defs::inds &inds) const;
    defs::inds flat_to_inds(const size_t &flat) const;
    size_t ndim_stride() const {return m_strides.size();}
    size_t get_dim(size_t idim) const {return m_shape[idim];}
    size_t get_last_dim() const {return get_dim(m_shape.size()-1);}
};


#endif //M7_INDEXER_H
