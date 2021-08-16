//
// Created by rja on 14/06/2020.
//

#include "Sparse.h"

void sparse::Network::resize(const size_t nrow) {
    DEBUG_ASSERT_GE(nrow, m_rows_icols.size(), "should be resizing to a larger number of rows");
    m_rows_icols.resize(nrow);
}

void sparse::Network::expand(const size_t delta_nrow) {
    resize(m_rows_icols.size() + delta_nrow);
}

size_t sparse::Network::nrow() const {
    return m_rows_icols.size();
}

void sparse::Network::add(const size_t &irow, const size_t &icol) {
    if (irow >= m_rows_icols.size()) {
        if (!m_resized_by_add) {
            log::warn("Resizing sparse matrix by adding a row (this entails reallocation which is inefficient)");
            log::warn("Call resize before filling if number of rows is known in advance");
            m_resized_by_add = true;
        }
        resize(irow + 1);
    }
    m_rows_icols[irow].push_front(icol);
}

bool sparse::Network::empty() { return m_rows_icols.empty(); }

const std::forward_list<size_t> &sparse::Network::operator[](const size_t &irow) {
    ASSERT(irow<nrow());
    return m_rows_icols[irow];
}
