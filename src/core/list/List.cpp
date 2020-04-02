//
// Created by Robert John Anderson on 2020-03-31.
//

#include "List.h"

List::List(size_t nsegment) : Table(nsegment), m_high_water_mark(nsegment, 0ul){}

const defs::inds &List::high_water_mark() const {
    return m_high_water_mark;
}

size_t List::push(const size_t &isegment) {
    assert(isegment<m_nsegment);
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark[isegment]++;
    if (tmp >= m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

size_t List::push(const size_t &isegment, const size_t &nrow) {
    assert(isegment<m_nsegment);
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark[isegment] += nrow;
    if (tmp >= m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

void List::zero() {
    // TODO: no need to memset zero here, only included initially for clarity in debugging
    Table::zero();
    m_high_water_mark.assign(0ul, m_nsegment);
}
