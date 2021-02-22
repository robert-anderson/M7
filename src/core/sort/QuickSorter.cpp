//
// Created by rja on 29/11/2020.
//

#include "QuickSorter.h"

Quicksorter::Quicksorter(Quicksorter::comp_t comp_fn) : m_comp_fn(comp_fn) {}

const size_t &Quicksorter::operator[](const size_t &i) const {
    ASSERT(i < m_inds.size());
    return m_inds[i];
}

void Quicksorter::sort(const size_t &hwm) {
    if (m_inds.size() < hwm) m_inds.reserve(hwm);
    m_inds.clear();
    // reset ordering
    for (size_t i = 0; i < hwm; ++i) m_inds.push_back(i);
    qs(0, hwm - 1);
    ASSERT(is_sorted(hwm));
}

void Quicksorter::sort(const TableBase &table) {
    sort(table.m_hwm);
}

bool Quicksorter::is_sorted(const size_t &hwm) {
    for (size_t irow = 1ul; irow < hwm; ++irow) {
        if (!m_comp_fn(m_inds[irow - 1], m_inds[irow])) return false;
    }
    return true;
}

bool Quicksorter::is_sorted(const TableBase &table) {
    return is_sorted(table.m_hwm);
}


void Quicksorter::swap(size_t ii1, size_t ii2) {
    auto i2 = m_inds[ii2];
    m_inds[ii2] = m_inds[ii1];
    m_inds[ii1] = i2;
}

size_t Quicksorter::partition(size_t iilo, size_t iihi) {
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

void Quicksorter::qs(size_t iilo, size_t iihi) {
    if (iihi != ~0ul && iilo < iihi) {
        auto iip = partition(iilo, iihi);
        qs(iilo, iip - 1);
        qs(iip + 1, iihi);
    }
}
