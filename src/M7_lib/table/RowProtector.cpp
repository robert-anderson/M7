//
// Created by rja on 19/05/2021.
//

#include "RowProtector.h"

RowProtector::RowProtector(const TableBase &table) : m_table(table), m_it(m_table.add_protector(this)) {
    on_resize(table.nrow());
}

RowProtector::~RowProtector() {
    m_table.erase_protector(this);
}

void RowProtector::protect(const size_t &irow) {
    // on_resize must be called when the source table is resized, otherwise we get OOB
    DEBUG_ASSERT_LT(irow, m_flags.size(), "row protection flag OOB");
    if (!m_flags[irow]) ++m_nprotected;
    m_flags[irow] = true;
}

void RowProtector::release(const size_t &irow) {
    DEBUG_ASSERT_LT(irow, m_flags.size(), "row protection flag OOB");
    if (m_flags[irow]) --m_nprotected;
    m_flags[irow] = false;
}

void RowProtector::on_resize(size_t nrow) {
    m_flags.resize(nrow);
}

bool RowProtector::is_protected(const size_t &irow) const {
    DEBUG_ASSERT_LT(irow, m_flags.size(), "row protection flag OOB");
    return m_flags[irow];
}

bool RowProtector::is_protected() const {
    return m_nprotected;
}