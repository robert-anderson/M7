//
// Created by rja on 29/11/2020.
//

#include "QuickSorter.h"

QuickSorter::QuickSorter(QuickSorter::comp_t comp_fn) : m_comp_fn(comp_fn) {}

const size_t &QuickSorter::operator[](const size_t &i) const {
    ASSERT(i < m_inds.size());
    return m_inds[i];
}

void QuickSorter::preserve_sort(const size_t &hwm) {
    if (m_inds.size() < hwm) m_inds.reserve(hwm);
    m_inds.clear();
    // reset ordering
    for (size_t i = 0; i < hwm; ++i) m_inds.push_back(i);
    qs(0, hwm - 1);
    ASSERT(is_preserve_sorted(hwm));
}

void QuickSorter::preserve_sort(const TableBase &table) {
    preserve_sort(table.m_hwm);
}

bool QuickSorter::is_preserve_sorted(const size_t &hwm) {
    for (size_t irow = 1ul; irow < hwm; ++irow) {
        if (!m_comp_fn(m_inds[irow - 1], m_inds[irow])) return false;
    }
    return true;
}

bool QuickSorter::is_preserve_sorted(const TableBase &table) {
    return is_preserve_sorted(table.m_hwm);
}


void QuickSorter::reorder_sort(TableBase &table) {
    qs(0, table.m_hwm - 1, table);
    ASSERT(is_reorder_sorted(table));
}

bool QuickSorter::is_reorder_sorted(const TableBase &table) {
    for (size_t irow = 1ul; irow < table.m_hwm; ++irow) {
        if (!m_comp_fn(irow - 1, irow)) return false;
    }
    return true;
}

void QuickSorter::swap(size_t ii1, size_t ii2) {
    auto i2 = m_inds[ii2];
    m_inds[ii2] = m_inds[ii1];
    m_inds[ii1] = i2;
}

size_t QuickSorter::partition(size_t iilo, size_t iihi) {
    auto ip = m_inds[iihi];
    auto ii = iilo - 1;

    for (size_t ij = iilo; ij <= iihi - 1; ij++) {
        if (m_comp_fn(m_inds[ij], ip)) {
            ii++;
            swap(ii, ij);
        }
    }
    swap(ii + 1, iihi);
    return ii + 1;
}

void QuickSorter::qs(size_t iilo, size_t iihi) {
    if (iihi != ~0ul && iilo < iihi) {
        auto iip = partition(iilo, iihi);
        qs(iilo, iip - 1);
        qs(iip + 1, iihi);
    }
}


void QuickSorter::swap(size_t ii1, size_t ii2, TableBase &table) {
    table.swap_rows(ii1, ii2);
}

size_t QuickSorter::partition(size_t iilo, size_t iihi, TableBase &table) {
    auto ip = iihi;
    auto ii = iilo - 1;

    for (size_t ij = iilo; ij <= iihi - 1; ij++) {
        if (m_comp_fn(ij, ip)) {
            ii++;
            swap(ii, ij, table);
        }
    }
    swap(ii + 1, iihi, table);
    return ii + 1;
}

void QuickSorter::qs(size_t iilo, size_t iihi, TableBase &table) {
    if (iihi != ~0ul && iilo < iihi) {
        auto iip = partition(iilo, iihi, table);
        qs(iilo, iip - 1, table);
        qs(iip + 1, iihi, table);
    }
}
