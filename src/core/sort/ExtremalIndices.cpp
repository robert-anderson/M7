//
// Created by rja on 29/11/2020.
//

#include <algorithm>
#include "ExtremalIndices.h"
#include "src/core/parallel/MPIAssert.h"

ExtremalIndices::ExtremalIndices(comparators::index_cmp_fn_t cmp_fn): m_cmp_fn(cmp_fn){}

const size_t &ExtremalIndices::nfound() const {
    return m_nfound;
}

size_t ExtremalIndices::nremain() const {
    DEBUG_ASSERT_LE(m_nfound, m_nind, "number found can't be more than the total number of elements!");
    return m_nind - m_nfound;
}

const size_t *ExtremalIndices::begin() const {
    return m_inds.data();
}

const size_t &ExtremalIndices::operator[](const size_t &ifound) const {
    ASSERT(ifound<m_nfound);
    return m_inds[ifound];
}

void ExtremalIndices::find(size_t nfind) {
    REQUIRE_NE(m_nind, ~0ul, "reset method must be called to initialise the number of indices to be sorted");
    auto limit = std::min(m_nind, m_nfound+nfind);
    std::partial_sort(m_inds.begin()+m_nfound, m_inds.begin()+limit, m_inds.end(), m_cmp_fn);
    m_nfound = limit;
}

void ExtremalIndices::reset(size_t hwm, defs::inds inds_ignore) {
    m_nfound = 0ul;
    m_nind = hwm - inds_ignore.size();
    if (m_inds.capacity() < m_nind) m_inds.reserve(m_nind);
    m_inds.clear();
    std::sort(inds_ignore.begin(), inds_ignore.end());
    // fill m_inds as an ordered array of indices. the indices are consecutive unless they appear in inds_ignore
    auto it_next_ignore=inds_ignore.begin();
    for (size_t i=0ul; i<m_nind; ++i){
        if (it_next_ignore!=inds_ignore.end() && i==*it_next_ignore) ++it_next_ignore;
        else m_inds.push_back(i);
    }
}

void ExtremalIndices::reset(const TableBase &table) {
    auto stack = table.m_free_rows;
    defs::inds irows_free;
    irows_free.reserve(stack.size());
    while(!stack.empty()) {
        irows_free.push_back(stack.top());
        stack.pop();
    }
    reset(table.nrow_nonzero(), irows_free);
}