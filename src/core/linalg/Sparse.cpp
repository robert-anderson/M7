//
// Created by rja on 14/06/2020.
//

#include "Sparse.h"

void sparse::Network::resize(const size_t& nrow) {
    if (nrow > m_rows_icols.size()) m_rows_icols.resize(nrow);
}

size_t sparse::Network::nrow() const {
    return m_rows_icols.size();
}

size_t sparse::Network::nentry() const {
    return m_nentry;
}

size_t sparse::Network::max_column_index() const {
    return m_max_icol;
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
    m_rows_icols[irow].push_back(icol);
    if (icol>m_max_icol) m_max_icol = icol;
    ++m_nentry;
}

void sparse::Network::checked_add(const size_t &irow, const size_t &icol) {
    const auto& row = m_rows_icols[irow];
    REQUIRE_FALSE(std::any_of(row.cbegin(), row.cend(), [&icol](const size_t& i){return i==icol;}),
                  "row entries in sparse network must have a unique column index");
    add(irow, icol);
}

void sparse::Network::add(const size_t &irow, const defs::inds &icols) {
    for (auto &icol: icols) add(irow, icol);
}

bool sparse::Network::empty() { return m_rows_icols.empty(); }

const defs::inds &sparse::Network::operator[](const size_t &irow) {
    ASSERT(irow<nrow());
    return m_rows_icols[irow];
}

std::vector<std::string> sparse::Network::row_to_strings(size_t irow) const {
    return utils::to_strings(m_rows_icols[irow]);
}

std::string sparse::Network::to_string() const {
    std::string out;
    for (size_t irow=0ul; irow<nrow(); ++irow) {
        if (m_rows_icols[irow].empty()) continue;
        out+= std::to_string(irow) + ": " + string_utils::join(row_to_strings(irow), ", ")+"\n";
    }
    return out;
}

sparse::Network sparse::Network::get_symmetrized() const {
    Network sym_net;
    sym_net.resize(nrow());
    REQUIRE_LT(max_column_index(), nrow(), "too many columns for this to be a symmetric matrix");
    for (size_t irow=0ul; irow<nrow(); ++irow) {
        const auto& icols = m_rows_icols[irow];
        for (size_t iicol=0ul; iicol < icols.size(); ++iicol) {
            const auto& icol = icols[iicol];
            sym_net.add(irow, icol);
            if (icol != irow) sym_net.checked_add(icol, irow);
        }
    }
    return sym_net;
}

void sparse::Network::get_row_subset(Network &subnet, size_t count, size_t displ) const {
    REQUIRE_LE(displ, nrow(), "row offset OOB");
    REQUIRE_LE(displ+count, nrow(), "row offset+count OOB");
    subnet.resize(count);
    auto begin = m_rows_icols.cbegin()+defs::inds::difference_type(displ);
    auto end = begin + defs::inds::difference_type(count);
    subnet.m_rows_icols = std::vector<defs::inds>(begin, end);
    // data is now copied, now update metadata
    for (const auto& row: subnet.m_rows_icols){
        if (row.empty()) continue;
        subnet.m_nentry+=row.size();
        auto max_element = std::max_element(row.cbegin(), row.cend());
        subnet.m_max_icol = std::max(subnet.m_max_icol, *max_element);
    }
}

sparse::Network sparse::Network::get_row_subset(size_t count, size_t displ) const {
    Network subnet;
    get_row_subset(subnet, count, displ);
    return subnet;
}
