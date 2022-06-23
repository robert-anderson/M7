//
// Created by Robert J. Anderson on 08/05/2021.
//

#include <numeric>
#include <M7_lib/io/Logging.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/sort/ExtremalIndices.h>

#include "TableBase.h"
#include "RowProtector.h"

TableBase::TableBase(size_t row_size) :
        m_bw(row_size), m_null_row_string(row_size, 0){}

TableBase::TableBase(const TableBase &other) :
        TableBase(other.m_bw.m_row_size){}

defs::buf_t *TableBase::begin() {
    return m_bw.m_begin;
}

const defs::buf_t *TableBase::begin() const {
    return m_bw.m_begin;
}

defs::buf_t *TableBase::begin(const size_t &irow) {
    return m_bw.m_begin + irow * row_size();
}

const defs::buf_t *TableBase::begin(const size_t &irow) const {
    return m_bw.m_begin + irow * row_size();
}

void TableBase::set_buffer(Buffer *buffer) {
    ASSERT(buffer);
    ASSERT(!m_bw.allocated())
    buffer->append_window(&m_bw);
}

bool TableBase::is_full() const {
    return m_hwm == nrow();
}

size_t TableBase::push_back(size_t nrow) {
    DEBUG_ASSERT_TRUE(row_size(), "cannot resize a table with zero row size");
    if (m_hwm + nrow > this->nrow()) expand(nrow);
    auto tmp = m_hwm;
    m_hwm += nrow;
    return tmp;
}

size_t TableBase::get_free_row() {
    if (m_free_rows.empty()) return push_back();
    auto irow = m_free_rows.top();
    m_free_rows.pop();
    return irow;
}

void TableBase::clear() {
    ASSERT(!is_protected());
    if (!m_bw.allocated()) return;
    std::memset(begin(), 0, row_size() * m_hwm);
    m_hwm = 0ul;
    while (!m_free_rows.empty()) m_free_rows.pop();
}

void TableBase::clear(const size_t &irow) {
    ASSERT(!is_protected(irow));
    std::memset(begin(irow), 0, row_size());
    m_free_rows.push(irow);
}

bool TableBase::is_cleared() const {
    return std::memcmp(begin(), m_null_row_string.data(), row_size()) == 0;
}

bool TableBase::is_cleared(const size_t &irow) const {
    return std::memcmp(begin(irow), m_null_row_string.data(), row_size()) == 0;
}

size_t TableBase::bw_size() const {
    return m_bw.m_size;
}

void TableBase::resize(size_t nrow, double factor) {
    DEBUG_ASSERT_TRUE(row_size(), "cannot resize, row size is zero");
    DEBUG_ASSERT_GE(nrow, m_hwm, "resize would discard uncleared data");
    m_bw.resize(nrow * row_size(), factor);
    DEBUG_ASSERT_LT(m_hwm, m_bw.m_size/row_size(), "resize has discarded uncleared data");
    for(const auto &rp : m_row_protectors) rp->on_resize(this->nrow());
}

void TableBase::expand(size_t nrow, double factor) {
    resize(this->nrow()+nrow, factor);
}

void TableBase::clear_rows(const defs::ivec_t &irows) {
    for (auto irow : irows) {
        clear(irow);
    }
}

void TableBase::post_insert(const size_t& iinsert) {}

void TableBase::insert_rows(const Buffer::Window &recv, size_t nrow, const std::list<recv_cb_t> &callbacks) {
    for (size_t irow_recv = 0; irow_recv < nrow; ++irow_recv) {
        auto irow = get_free_row();
        std::memcpy(begin(irow), recv.m_begin + irow_recv * row_size(), row_size());
        post_insert(irow);
        for (auto f: callbacks) f(irow);
    }
}
void TableBase::transfer_rows(const defs::ivec_t &irows, size_t irank_send, size_t irank_recv, const std::list<recv_cb_t>& callbacks){
    DEBUG_ASSERT_NE_ALL(irank_recv, irank_send, "sending and recving ranks should never be the same");
    if (!m_transfer) m_transfer = smart_ptr::make_unique<RowTransfer>(m_bw.name());
    size_t nrow = 0;
    if (mpi::i_am(irank_send)){
        auto& send_bw = m_transfer->m_send_bw;
        nrow = irows.size();
        mpi::send(&nrow, 1, irank_recv, m_transfer->m_nrow_p2p_tag);
        if (!nrow){
            log::debug_("Sending rank notifying recving rank that no rows are transferred");
            return;
        }
        log::info_("Transferring {} rows outward to rank {}", nrow, irank_recv);

        auto size_required = nrow*row_size();
        if (send_bw.m_size < size_required) send_bw.resize(size_required);
        for (auto iirow = 0ul; iirow < nrow; ++iirow) {
            const auto &irow = irows[iirow];
            std::memcpy(send_bw.m_begin + iirow * row_size(), begin(irow), row_size());
        }
        mpi::send(send_bw.m_begin, row_size() * nrow, irank_recv, m_transfer->m_irows_p2p_tag);
        /*
         * sent rows can now be erased
         */
        clear_rows(irows);
    }
    if (mpi::i_am(irank_recv)){
        auto& recv_bw = m_transfer->m_send_bw;
        mpi::recv(&nrow, 1, irank_send, m_transfer->m_nrow_p2p_tag);
        if (!nrow){
            log::debug_("Recving rank notified by sending rank that no rows are transferred");
            return;
        }
        auto size_required = nrow * row_size();
        if (recv_bw.m_size<size_required) recv_bw.resize(size_required);
        log::info_("Transferring {} rows inward from rank {}", nrow, irank_send);
        mpi::recv(recv_bw.m_begin, row_size() * nrow, irank_send, m_transfer->m_irows_p2p_tag);
        /*
         * now emplace received rows in TableBase buffer window, and call all callbacks for each
         */
        insert_rows(recv_bw, nrow, callbacks);
    }
}

void TableBase::copy_row_in(const TableBase &src, size_t irow_src, size_t irow_dst) {
    ASSERT(irow_dst < m_hwm);
    std::memcpy(begin(irow_dst), src.begin(irow_src), row_size());
}

void TableBase::swap_rows(const size_t &irow, const size_t &jrow) {
    if (irow == jrow) return;
    auto iptr = begin(irow);
    auto jptr = begin(jrow);
    std::swap_ranges(iptr, iptr+row_size(), jptr);
}

std::string TableBase::to_string(const defs::ivec_t *ordering) const {
    std::string out;
    auto row_ptr = begin();
    for (size_t irow=0ul; irow<m_hwm; ++irow){
        for (size_t idword=0ul; idword<row_size(); ++idword){
            out+=std::to_string(static_cast<int>(row_ptr[idword]))+" ";
        }
        row_ptr+=row_size();
        out+="\n";
    }
    return out;
}

void TableBase::all_gatherv(const TableBase &src) {
    clear();
    defs::ivec_t nrows(mpi::nrank());
    defs::ivec_t counts(mpi::nrank());
    defs::ivec_t displs(mpi::nrank());
    DEBUG_ASSERT_EQ(src.row_size(), row_size(),
                    "the size of rows being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.m_hwm, nrows);
    counts = nrows;
    for (auto &v: counts) v *= row_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrow_total = std::accumulate(nrows.cbegin(), nrows.cend(), 0ul);
    push_back(nrow_total);
    mpi::all_gatherv(src.begin(), src.m_hwm * row_size(), begin(), counts, displs);
    post_insert_range(0, nrow_total);
}

void TableBase::gatherv(const TableBase &src, size_t irank) {
    if (mpi::i_am(irank)) clear();
    defs::ivec_t nrows(mpi::nrank());
    defs::ivec_t counts(mpi::nrank());
    defs::ivec_t displs(mpi::nrank());
    DEBUG_ASSERT_EQ(src.row_size(), row_size(),
                    "the size of rows being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.m_hwm, nrows);
    counts = nrows;
    for (auto &v: counts) v *= row_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrow_total = std::accumulate(nrows.cbegin(), nrows.cend(), 0ul);
    if (mpi::i_am(irank))
        push_back(nrow_total);
    mpi::gatherv(src.begin(), src.m_hwm * row_size(), begin(), counts, displs, irank);
    if (mpi::i_am(irank))
        post_insert_range(0, nrow_total);
}

void TableBase::erase_protector(RowProtector *rp) const {
    m_row_protectors.erase(rp->m_it);
}

bool TableBase::is_protected() const {
    return std::any_of(m_row_protectors.cbegin(), m_row_protectors.cend(),
                       [](const RowProtector* rp){return rp->is_protected();});
}

bool TableBase::is_protected(const size_t& irow) const {
    return std::any_of(m_row_protectors.cbegin(), m_row_protectors.cend(),
                       [&](const RowProtector* rp){return rp->is_protected(irow);});
}

size_t TableBase::nrow_nonzero() const {
    return m_hwm-m_free_rows.size();
}

TableBase::Loc::Loc(size_t irank, size_t irow) : m_irank(irank), m_irow(irow){
#ifndef NDEBUG
    mpi::bcast(irank);
    mpi::bcast(irow);
    DEBUG_ASSERT_EQ(m_irank, irank, "rank index in TableBase::Loc should be consistent across all ranks");
    DEBUG_ASSERT_EQ(m_irow, irow, "row index in TableBase::Loc should be consistent across all ranks");
#endif
}

TableBase::Loc::operator bool() const {
    return m_irank!=~0ul;
}

bool TableBase::Loc::is_mine() const {
    return mpi::i_am(m_irank);
}

bool TableBase::Loc::operator==(const TableBase::Loc &other) {
    return m_irank==other.m_irank and m_irow==other.m_irow;
}

bool TableBase::Loc::operator!=(const TableBase::Loc &other) {
    return !(*this==other);
}
