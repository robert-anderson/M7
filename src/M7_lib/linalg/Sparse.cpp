//
// Created by Robert J. Anderson on 14/06/2020.
//

#include "Sparse.h"
#include "M7_lib/util/String.h"

void sparse::dynamic::Network::resize(uint_t nrow) {
    if (nrow > m_rows_icols.size()) m_rows_icols.resize(nrow);
}

uint_t sparse::dynamic::Network::nrow() const {
    return m_rows_icols.size();
}

uint_t sparse::dynamic::Network::nentry() const {
    return m_nentry;
}

uint_t sparse::dynamic::Network::nentry(uint_t irow) const {
    return (*this)[irow].size();
}

uint_t sparse::dynamic::Network::max_column_index() const {
    return m_max_icol;
}

uint_t sparse::dynamic::Network::add(uint_t irow, uint_t icol) {
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
    return m_rows_icols[irow].size()-1;
}

uint_t sparse::dynamic::Network::insert(uint_t irow, uint_t icol) {
    auto& row = m_rows_icols[irow];
    auto element = std::find(row.begin(), row.end(), icol);
    if (element==row.end()) return add(irow, icol);
    *element = icol;
    return std::distance(row.begin(), element);
}

void sparse::dynamic::Network::add(uint_t irow, const uintv_t &icols) {
    for (auto &icol: icols) add(irow, icol);
}

void sparse::dynamic::Network::insert(uint_t irow, const uintv_t &icols) {
    for (auto &icol: icols) insert(irow, icol);
}

bool sparse::dynamic::Network::empty() const {
    return m_rows_icols.empty();
}

bool sparse::dynamic::Network::empty(uint_t irow) const {
    return (*this)[irow].empty();
}

const uintv_t &sparse::dynamic::Network::operator[](uint_t irow) const{
    ASSERT(irow<nrow());
    return m_rows_icols[irow];
}

str_t sparse::dynamic::Network::row_to_string(uint_t irow) const {
    return convert::to_string(m_rows_icols[irow]);
}

str_t sparse::dynamic::Network::to_string() const {
    str_t out;
    for (uint_t irow=0ul; irow<nrow(); ++irow) {
        if (m_rows_icols[irow].empty()) continue;
        out+= std::to_string(irow) + ": " + row_to_string(irow)+"\n";
    }
    return out;
}

sparse::dynamic::Network sparse::dynamic::Network::get_symmetrized() const {
    Network sym_net;
    sym_net.resize(nrow());
    REQUIRE_LT(max_column_index(), nrow(), "too many columns for this to be a symmetric matrix");
    for (uint_t irow=0ul; irow<nrow(); ++irow) {
        const auto& icols = m_rows_icols[irow];
        for (uint_t iicol=0ul; iicol < icols.size(); ++iicol) {
            const auto& icol = icols[iicol];
            sym_net.add(irow, icol);
            if (icol != irow) sym_net.insert(icol, irow);
        }
    }
    return sym_net;
}

void sparse::dynamic::Network::get_row_subset(Network &subnet, uint_t count, uint_t displ) const {
    REQUIRE_LE(displ, nrow(), "row offset OOB");
    REQUIRE_LE(displ+count, nrow(), "row offset+count OOB");
    subnet.resize(count);
    auto begin = m_rows_icols.cbegin()+uintv_t::difference_type(displ);
    auto end = begin + uintv_t::difference_type(count);
    subnet.m_rows_icols = v_t<uintv_t>(begin, end);
    // data is now copied, now update metadata
    for (const auto& row: subnet.m_rows_icols){
        if (row.empty()) continue;
        subnet.m_nentry+=row.size();
        auto max_element = std::max_element(row.cbegin(), row.cend());
        subnet.m_max_icol = std::max(subnet.m_max_icol, *max_element);
    }
}

sparse::dynamic::Network sparse::dynamic::Network::get_row_subset(uint_t count, uint_t displ) const {
    Network subnet;
    get_row_subset(subnet, count, displ);
    return subnet;
}

uintv_t sparse::fixed::Base::make_counts(const sparse::dynamic::Network &src) {
    uintv_t counts;
    counts.reserve(src.nrow());
    for (uint_t irow=0ul; irow<src.nrow(); ++irow) counts.push_back(src.nentry(irow));
    return counts;
}

uintv_t sparse::fixed::Base::make_displs(const uintv_t &counts) {
    uintv_t displs;
    /*
     * include the last one so that the cend iterator for the last entry is directly accessible
     */
    displs.reserve(counts.size()+1);
    displs.push_back(0ul);
    for (auto it = counts.cbegin(); it!=counts.cend(); ++it) displs.push_back(displs.back()+*it);
    return displs;
}

sparse::fixed::Base::Base(const uintv_t &counts) :
    m_nrow(counts.size()), m_max_nentry(counts.empty() ? 0ul : *std::max(counts.cbegin(), counts.cend())),
    m_displs(make_displs(counts)), m_nentry(m_displs.back()) {}

sparse::fixed::Base::Base(const sparse::dynamic::Network &src) : Base(make_counts(src)){
    DEBUG_ASSERT_EQ(m_nentry, src.nentry(),
                    "number of entries calculated from offsets should match that counted by the dynamic::Network instance");
    DEBUG_ASSERT_EQ(m_max_nentry, src.max_column_index()+1,
                    "max number of entries in row calculated from offsets should match dynamic::Network::m_icol_max+1");
}

uintv_t sparse::fixed::Network::make_entries(const sparse::dynamic::Network &src) {
    uintv_t out;
    out.reserve(src.nrow());
    for (uint_t irow=0ul; irow<src.nrow(); ++irow)
        out.insert(out.cend(), src[irow].cbegin(), src[irow].cend());
    return out;
}

uintv_t::const_iterator sparse::fixed::Network::citer(uint_t irow) const {
    auto it = m_entries.cbegin();
    std::advance(it, m_displs[irow]);
    return it;
}

sparse::fixed::Network::Network(const sparse::dynamic::Network &src) : Base(src), m_entries(make_entries(src)){
    DEBUG_ASSERT_EQ(m_entries.size(), m_nentry, "incorrect number of entries");
}

uintv_t::const_iterator sparse::fixed::Network::cbegin(uint_t irow) const {
    if (irow>=m_nrow) return m_entries.cend();
    return citer(irow);
}

uintv_t::const_iterator sparse::fixed::Network::cend(uint_t irow) const {
    if (irow>=m_nrow) return m_entries.cend();
    return citer(irow+1);
}

uint_t sparse::fixed::Network::nentry(uint_t irow) const {
    return std::distance(cbegin(irow), cend(irow));
}

sparse::inverse::Network::Network(const sparse::dynamic::Network &src) {
    for (auto irow=0ul; irow<src.nrow(); ++irow)
        for (const auto& it: src[irow])
            m_conns.insert({irow, it});
}

sparse::inverse::Network::Network(const sparse::fixed::Network &src) {
    for (auto irow=0ul; irow<src.m_nrow; ++irow)
        for (auto it = src.cbegin(irow); it != src.cend(irow); ++it) m_conns.insert({irow, *it});
}

bool sparse::inverse::Network::lookup(uint_t i, uint_t j) const {
    return m_conns.find({i, j}) != m_conns.cend();
}
