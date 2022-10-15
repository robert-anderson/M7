//
// Created by Robert J. Anderson on 29/11/2020.
//

#include "LambdaQuickSorter.h"

LambdaQuickSorter::LambdaQuickSorter(comparators::index_cmp_fn_t cmp_fn) : m_cmp_fn(cmp_fn) {}

const uint_t &LambdaQuickSorter::operator[](const uint_t &i) const {
    ASSERT(i < m_inds.size());
    return m_inds[i];
}

void LambdaQuickSorter::preserve_sort(const uint_t &hwm) {
    if (m_inds.size() < hwm) m_inds.reserve(hwm);
    m_inds.clear();
    // reset ordering
    for (uint_t i = 0; i < hwm; ++i) m_inds.push_back(i);
    qs(0, hwm - 1);
    ASSERT(is_preserve_sorted(hwm));
}

void LambdaQuickSorter::preserve_sort(const TableBase &table) {
    preserve_sort(table.nrow_in_use());
}

bool LambdaQuickSorter::is_preserve_sorted(const uint_t &hwm) {
    for (uint_t irow = 1ul; irow < hwm; ++irow) {
        if (m_cmp_fn(m_inds[irow], m_inds[irow-1]) &&!m_cmp_fn(m_inds[irow - 1], m_inds[irow])) return false;
    }
    return true;
}

bool LambdaQuickSorter::is_preserve_sorted(const TableBase &table) {
    return is_preserve_sorted(table.nrow_in_use());
}


void LambdaQuickSorter::reorder_sort(TableBase &table) {
    qs(0, table.nrow_in_use() - 1, table);
    ASSERT(is_reorder_sorted(table));
}

bool LambdaQuickSorter::is_reorder_sorted(const TableBase &table) {
    for (uint_t irow = 1ul; irow < table.nrow_in_use(); ++irow) {
        if (m_cmp_fn(irow, irow-1) && !m_cmp_fn(irow - 1, irow)) return false;
    }
    return true;
}

void LambdaQuickSorter::swap(uint_t ii1, uint_t ii2) {
    auto i2 = m_inds[ii2];
    m_inds[ii2] = m_inds[ii1];
    m_inds[ii1] = i2;
}

uint_t LambdaQuickSorter::partition(uint_t iilo, uint_t iihi) {
    auto ip = m_inds[iihi];
    auto ii = iilo - 1;

    for (uint_t ij = iilo; ij <= iihi - 1; ij++) {
        if (m_cmp_fn(m_inds[ij], ip)) {
            ii++;
            swap(ii, ij);
        }
    }
    swap(ii + 1, iihi);
    return ii + 1;
}

void LambdaQuickSorter::qs(uint_t iilo, uint_t iihi) {
    if (iihi != ~0ul && iilo < iihi) {
        auto iip = partition(iilo, iihi);
        qs(iilo, iip - 1);
        qs(iip + 1, iihi);
    }
}


void LambdaQuickSorter::swap(uint_t ii1, uint_t ii2, TableBase &table) {
    table.swap_records(ii1, ii2);
}

uint_t LambdaQuickSorter::partition(uint_t iilo, uint_t iihi, TableBase &table) {
    auto ip = iihi;
    auto ii = iilo - 1;

    for (uint_t ij = iilo; ij <= iihi - 1; ij++) {
        if (m_cmp_fn(ij, ip)) {
            ii++;
            swap(ii, ij, table);
        }
    }
    swap(ii + 1, iihi, table);
    return ii + 1;
}

void LambdaQuickSorter::qs(uint_t iilo, uint_t iihi, TableBase &table) {
    if (iihi != ~0ul && iilo < iihi) {
        auto iip = partition(iilo, iihi, table);
        qs(iilo, iip - 1, table);
        qs(iip + 1, iihi, table);
    }
}
