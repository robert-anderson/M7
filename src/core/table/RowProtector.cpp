//
// Created by rja on 19/05/2021.
//

#include "RowProtector.h"

void RowProtector::protect(const size_t &irow) {
    if (!m_flags[irow]) ++m_nprotected;
    m_flags[irow] = true;
}

void RowProtector::release(const size_t &irow) {
    if (m_flags[irow]) --m_nprotected;
    m_flags[irow] = false;
}

void RowProtector::on_resize(size_t nrow) {
    m_flags.resize(nrow);
}

bool RowProtector::is_protected(const size_t &irow) const {
    ASSERT(irow < m_flags.size());
    return m_flags[irow];
}

bool RowProtector::is_protected() const {
    return m_nprotected;
}