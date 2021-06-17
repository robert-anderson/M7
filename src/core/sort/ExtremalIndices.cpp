//
// Created by rja on 29/11/2020.
//

#include <algorithm>
#include "ExtremalIndices.h"
#include "src/core/parallel/MPIAssert.h"

ExtremalIndices::ExtremalIndices(ExtremalIndices::cmp_fn_t cmp_fn) :
    // reverse arg order to match ordering of QuickSort implementation
    m_cmp_fn([cmp_fn](const size_t& i2, const size_t& i1){return cmp_fn(i1, i2);}){}

const size_t &ExtremalIndices::nfound() const {
    return m_nfound;
}

const size_t *ExtremalIndices::begin() const {
    return m_inds.data()+(m_inds.size()-nfound());
}

const size_t &ExtremalIndices::operator[](const size_t &ifound) const {
    ASSERT(ifound<m_nfound);
    /*
     * the results should be read from the back of the heap vector
     */
    return m_inds[m_inds.size()-1-ifound];
}

void ExtremalIndices::find(size_t nfind) {
    REQUIRE_NE(m_hwm, ~0ul, "reset method must be called to initialise high water mark");
    nfind+=m_nfound;
    for (; m_nfound<std::min(m_hwm, nfind); ++m_nfound)
        std::pop_heap(m_inds.begin(), m_inds.end()-m_nfound, m_cmp_fn);
}

void ExtremalIndices::reset(size_t hwm) {
    m_nfound = 0ul;
    m_hwm = hwm;
    if (m_inds.capacity() < m_hwm) m_inds.reserve(m_hwm);
    m_inds.clear();
    for (size_t i = 0ul; i < m_hwm; ++i) m_inds.push_back(i);
    std::make_heap(m_inds.begin(), m_inds.end(), m_cmp_fn);
}

void ExtremalIndices::reset(const TableBase &table) { reset(table.m_hwm); }