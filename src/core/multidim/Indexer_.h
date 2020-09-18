//
// Created by RJA on 18/09/2020.
//

#ifndef M7_INDEXER__H
#define M7_INDEXER__H

#include "ArrayFormatBase.h"


template<size_t nind, size_t depth>
struct IndexerBase {
    const size_t m_offset;
    const ArrayFormatBase<nind> &m_format;

    IndexerBase(const ArrayFormatBase<nind> &format, const size_t &i) :
            m_offset(format.strides()[0]*i), m_format(format) {
        static_assert(depth == 1, "This ctor should only be called by the ArrayFormat [] operator");
    }

    IndexerBase(const IndexerBase<nind, depth - 1> &parent, const size_t &i) :
            m_offset(parent.m_offset + parent.m_format.strides()[depth-1] * i), m_format(parent.m_format) {
        static_assert(depth > 0ul, "depth parameter arg should always be nonzero");
        assert(i < parent.m_format.shape()[depth-1]);
    }

    size_t nelement() const {
        return m_format.strides()[depth-1] * m_format.shape()[depth-1];
    }

    const size_t& offset() const {
        return m_offset;
    }
};

template<size_t nind, size_t depth>
struct Indexer_ : IndexerBase<nind, depth> {
    Indexer_(const ArrayFormatBase<nind> &format, const size_t &i) : IndexerBase<nind, depth>(format, i) {}

    Indexer_(const Indexer_<nind, depth - 1> &parent, const size_t &i) : IndexerBase<nind, depth>(parent, i) {}

    Indexer_<nind, depth + 1> operator[](size_t i) {
        static_assert(depth <= nind, "depth of indexing may not exceed number of indices");
        return Indexer_<nind, depth + 1>(*this, i);
    }

};

template<size_t nind>
struct Indexer_<nind, nind> : IndexerBase<nind, nind> {
    Indexer_(const ArrayFormatBase<nind> &format, const size_t &i) : IndexerBase<nind, 1ul>(format, i) {}

    Indexer_(const IndexerBase<nind, nind - 1> &parent, const size_t &i) : IndexerBase<nind, nind>(parent, i) {}

    operator size_t() const { return IndexerBase<nind, nind>::m_offset; }
};



#endif //M7_INDEXER__H
