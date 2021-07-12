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
    DEBUG_ASSERT_LE(m_nfound, m_hwm, "number found can't be more than the total number of elements!");
    return m_hwm - m_nfound;
}

const size_t *ExtremalIndices::begin() const {
    return m_inds.data();
}

const size_t &ExtremalIndices::operator[](const size_t &ifound) const {
    ASSERT(ifound<m_nfound);
    return m_inds[ifound];
}

void ExtremalIndices::find(size_t nfind) {
    REQUIRE_NE(m_hwm, ~0ul, "reset method must be called to initialise high water mark");
    auto limit = std::min(m_hwm, nfind+m_nfound);
    std::partial_sort(m_inds.begin()+m_nfound, m_inds.begin()+limit, m_inds.end(), m_cmp_fn);
    m_nfound = limit;
}

void ExtremalIndices::reset(size_t hwm) {
    m_nfound = 0ul;
    m_hwm = hwm;
    if (m_inds.capacity() < m_hwm) m_inds.reserve(m_hwm);
    m_inds.clear();
    // fill m_inds as an ordered, consecutive array of indices
    std::iota(m_inds.begin(), m_inds.end(), 0ul);
}

void ExtremalIndices::reset(const TableBase &table) { reset(table.m_hwm); }