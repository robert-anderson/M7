//
// Created by rja on 29/11/2020.
//

#include <algorithm>
#include "ExtremalValues.h"
#include "src/core/parallel/MPIAssert.h"

ExtremalValues::ExtremalValues(ExtremalValues::comp_t comp_fn) :
    // reverse arg order to match ordering of QuickSort implementation
    m_comp_fn([comp_fn](const size_t& i2, const size_t& i1){return comp_fn(i1, i2);}){}

const size_t &ExtremalValues::nfound() const {
    return m_nfound;
}

const size_t *ExtremalValues::begin() const {
    return m_inds.data()+(m_inds.size()-nfound());
}

const size_t &ExtremalValues::operator[](const size_t &ifound) const {
    ASSERT(ifound<m_nfound);
    /*
     * the results should be read from the back of the heap vector
     */
    return m_inds[m_inds.size()-1-ifound];
}

void ExtremalValues::find(size_t nfind) {
    MPI_REQUIRE(m_hwm!=~0ul, "reset method must be called to initialise high water mark");
    nfind+=m_nfound;
    for (; m_nfound<std::min(m_hwm, nfind); ++m_nfound)
        std::pop_heap(m_inds.begin(), m_inds.end()-m_nfound, m_comp_fn);
}

void ExtremalValues::reset(size_t hwm) {
    m_nfound = 0ul;
    m_hwm = hwm;
    if (m_inds.capacity() < m_hwm) m_inds.reserve(m_hwm);
    m_inds.clear();
    for (size_t i = 0ul; i < m_hwm; ++i) m_inds.push_back(i);
    std::make_heap(m_inds.begin(), m_inds.end(), m_comp_fn);
}

void ExtremalValues::reset(const TableBase &table) { reset(table.m_hwm); }

const defs::inds &ExtremalValues::inds() const {
    return m_inds;
}
