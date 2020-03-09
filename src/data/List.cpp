//
// Created by Robert John Anderson on 2020-02-20.
//

#include "List.h"

List::List(defs::data_t *data_external) : Table(data_external) {}

size_t List::high_water_mark() const {
    return m_high_water_mark;
}

size_t List::push() {
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark++;
    if (tmp>=m_nrow) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

size_t List::push(const size_t &nrow) {
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark += nrow;
    if (tmp>=m_nrow) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

void List::zero() {
    // TODO: no need to memset zero here, only included initially for clarity in debugging
    Table::zero();
    m_high_water_mark = 0ul;
}

void List::print(size_t irank) const {
    Table::print(high_water_mark(), irank);
}

void List::high_water_mark(const size_t &value) {
    m_high_water_mark = value;
}
