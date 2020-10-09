//
// Created by RJA on 18/09/2020.
//

#ifndef M7_INDEXER_H
#define M7_INDEXER_H

#include "NdSpecifierBase.h"


template<typename selector_t, size_t nind, size_t nind_unspec>
struct IndexerBase {

    // the index of the shape array element being specermined in this object
    static constexpr auto nind_spec = nind-nind_unspec;

    const size_t m_offset;
    NdSpecifierBase<selector_t, nind> &m_spec;

    /*
     * TODO: move this to special case partial template where nind_unspec == nind
     */
    IndexerBase(NdSpecifierBase<selector_t, nind> &spec, const size_t &i) :
            m_offset(spec.strides()[0]*i), m_spec(spec) {
        static_assert(nind_spec == 1ul, "This ctor should only be called by the NdSpecifier [] operator");
    }

    IndexerBase(const IndexerBase<selector_t, nind, nind_unspec + 1> &parent, const size_t &i) :
            m_offset(parent.m_offset + parent.m_spec.strides()[nind_spec] * i), m_spec(parent.m_spec) {
        static_assert(nind_spec > 0ul, "number of specermined indices should always be nonzero when constructing from parent IndexerBase");
        assert(i < parent.m_spec.shape()[nind_spec]);
    }

    size_t nelement() const {
        return m_spec.strides()[nind_spec] * m_spec.shape()[nind_spec];
    }

    const size_t& offset() const {
        return m_offset;
    }
};

template<typename selector_t, size_t nind, size_t nind_unspec>
struct Indexer : IndexerBase<selector_t, nind, nind_unspec> {
    Indexer(NdSpecifierBase<selector_t, nind> &spec, const size_t &i) :
    IndexerBase<selector_t, nind, nind_unspec>(spec, i) {}

    Indexer(const Indexer<selector_t, nind, nind_unspec> &parent, const size_t &i) :
    IndexerBase<selector_t, nind, nind_unspec>(parent, i) {}

    using IndexerBase<selector_t, nind, nind_unspec>::nind_spec;
    Indexer<selector_t, nind, nind_unspec - 1> operator[](const size_t& i) {
        static_assert(nind_spec < nind, "depth of indexing may not exceed number of indices");
        return Indexer<selector_t, nind, nind_unspec - 1>(*this, i);
    }

};

template<typename selector_t, size_t nind>
struct Indexer<selector_t, nind, 1ul> : IndexerBase<selector_t, nind, 1ul> {

    typedef IndexerBase<selector_t, nind, 1ul> base_t;
    Indexer(NdSpecifierBase<selector_t, nind> &spec, const size_t &i) :
            base_t(spec, i) {}

    Indexer(const Indexer<selector_t, nind, 2ul> &parent, const size_t &i) :
            base_t(parent, i) {}

    using base_t::nind_spec;
            /*
    Indexer<selector_t, nind, depth + 1> operator[](const size_t& i) {
        static_assert(depth <= nind, "depth of indexing may not exceed number of indices");
        return Indexer<selector_t, nind, depth + 1>(*this, i);
    }*/

    using base_t::m_offset;
    using base_t::m_spec;
    typename selector_t::accessor_t operator[](const size_t& i) {
        return m_spec.select(m_offset+i);
    }


};

//template<typename selector_t, size_t nind>
//struct Indexer<selector_t, nind, nind> : IndexerBase<selector_t, nind, nind> {
//    Indexer(const ArrayspecBase<selector_t, nind> &spec, const size_t &i) :
//    IndexerBase<selector_t, nind, 1ul>(spec, i) {}
//
//    Indexer(const IndexerBase<selector_t, nind, nind - 1> &parent, const size_t &i) :
//    IndexerBase<selector_t, nind, nind>(parent, i) {}
//
//    //operator size_t() const { return IndexerBase<nind, nind>::m_offset; }
//
//    using IndexerBase<selector_t, nind, nind>::m_offset;
//    using IndexerBase<selector_t, nind, nind>::m_spec;
//
//    typename selector_t::accessor_t select(){
//        return m_spec.select(m_offset);
//    }
//    /*
//    typename selector_t::const_accessor_t select() const{
//        return m_spec.select(m_offset);
//    }
//     */
//};



#endif //M7_INDEXER_H
