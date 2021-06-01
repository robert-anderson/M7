//
// Created by rja on 08/05/2021.
//

#include "TableBase.h"
#include "src/core/io/Logging.h"
#include "src/core/parallel/MPIAssert.h"
#include "src/core/sort/ExtremalIndices.h"
#include "RowProtector.h"

TableBase::TableBase(size_t row_dsize) :
        m_row_dsize(row_dsize), m_row_size(row_dsize * defs::nbyte_data){
    // row must have non-zero size
    ASSERT(row_dsize);
}

TableBase::TableBase(const TableBase &other) :
        TableBase(other.m_row_dsize){}

defs::data_t *TableBase::dbegin() {
    return m_bw.m_dbegin;
}

const defs::data_t *TableBase::dbegin() const {
    return m_bw.m_dbegin;
}

defs::data_t *TableBase::dbegin(const size_t &irow) {
    return m_bw.m_dbegin + irow * m_row_dsize;
}

const defs::data_t *TableBase::dbegin(const size_t &irow) const {
    return m_bw.m_dbegin + irow * m_row_dsize;
}

void TableBase::set_buffer(Buffer *buffer) {
    ASSERT(buffer);
    ASSERT(!m_bw.allocated())
    buffer->append_window(&m_bw);
}

bool TableBase::is_full() const {
    return m_hwm == m_nrow;
}

size_t TableBase::push_back(size_t nrow) {
    if (m_hwm + nrow > m_nrow) expand(nrow);
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
    std::memset(dbegin(), 0, m_row_size * m_hwm);
    m_hwm = 0ul;
    while (!m_free_rows.empty()) m_free_rows.pop();
}

void TableBase::clear(const size_t &irow) {
    ASSERT(!is_protected(irow));
    std::memset(dbegin(irow), 0, m_row_size);
    m_free_rows.push(irow);
}

bool TableBase::is_cleared() const {
    return std::all_of(dbegin(), dbegin() + m_row_dsize * m_nrow, [](const defs::data_t &i) { return i == 0; });
}

bool TableBase::is_cleared(const size_t &irow) const {
    return std::all_of(dbegin(irow), dbegin(irow) + m_row_dsize, [](const defs::data_t &i) { return i == 0; });
}

size_t TableBase::bw_dsize() const {
    return m_bw.dsize();
}

void TableBase::resize(size_t nrow) {
    assert(nrow > m_nrow);
    m_bw.resize(nrow * m_row_dsize);
    for(const auto &rp : m_row_protectors) rp->on_resize(nrow);
    m_nrow = nrow;
}

void TableBase::expand(size_t nrow) {
    resize(m_nrow + nrow);
}

void TableBase::expand(size_t nrow, double expansion_factor) {
    resize((m_nrow + nrow)*(1+expansion_factor));
}

void TableBase::erase_rows(const defs::inds &irows) {
    for (auto irow : irows) {
        clear(irow);
    }
}

void TableBase::post_insert(const size_t& iinsert) {}

void TableBase::insert_rows(const Buffer::Window &recv, size_t nrow, const std::list<recv_cb_t> &callbacks) {
    for (size_t irow_recv = 0; irow_recv < nrow; ++irow_recv) {
        auto irow = get_free_row();
        std::memcpy(dbegin(irow), recv.m_dbegin + irow_recv * m_row_dsize, m_row_size);
        post_insert(irow);
        for (auto f: callbacks) f(irow);
    }
}
void TableBase::transfer_rows(const defs::inds &irows, size_t irank_send, size_t irank_recv, const std::list<recv_cb_t>& callbacks){
    MPI_ASSERT_ALL(irank_recv!=irank_send, "sending and recving ranks should never be the same");
    if (!m_transfer) m_transfer = std::unique_ptr<RowTransfer>(new RowTransfer(m_bw.name()));
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

        send_bw.make_room(nrow * m_row_dsize);
        for (auto iirow = 0ul; iirow < nrow; ++iirow) {
            const auto &irow = irows[iirow];
            std::memcpy(send_bw.m_dbegin + iirow * m_row_dsize, dbegin(irow), m_row_size);
        }
        mpi::send(send_bw.m_dbegin, m_row_dsize * nrow, irank_recv, m_transfer->m_irows_p2p_tag);
        /*
         * sent rows can now be erased
         */
        erase_rows(irows);
    }
    if (mpi::i_am(irank_recv)){
        auto& recv_bw = m_transfer->m_send_bw;
        mpi::recv(&nrow, 1, irank_send, m_transfer->m_nrow_p2p_tag);
        if (!nrow){
            log::debug_("Recving rank notified by sending rank that no rows are transferred");
            return;
        }
        recv_bw.make_room(nrow * m_row_dsize);
        log::info_("Transferring {} rows inward from rank {}", nrow, irank_send);
        mpi::recv(recv_bw.m_dbegin, m_row_dsize * nrow, irank_send, m_transfer->m_irows_p2p_tag);
        /*
         * now emplace received rows in TableBase buffer window, and call all callbacks for each
         */
        insert_rows(recv_bw, nrow, callbacks);
    }
}

void TableBase::copy_row_in(const TableBase &src, size_t irow_src, size_t irow_dst) {
    ASSERT(irow_dst < m_hwm);
    std::memcpy(dbegin(irow_dst), src.dbegin(irow_src), m_row_size);
}

void TableBase::swap_rows(const size_t &irow, const size_t &jrow) {
    if (irow == jrow) return;
    auto iptr = dbegin(irow);
    auto jptr = dbegin(jrow);
    for (size_t idword = 0ul; idword < m_row_dsize; ++idword) {
        std::swap(*iptr, *jptr);
        ++iptr;
        ++jptr;
    }
}

std::string TableBase::to_string(const defs::inds *ordering) const {
    return "";
}

void TableBase::all_gatherv(const TableBase &src) {
    clear();
    defs::inds nrows(mpi::nrank());
    defs::inds counts(mpi::nrank());
    defs::inds displs(mpi::nrank());
    ASSERT(src.m_row_dsize == m_row_dsize);
    mpi::all_gather(src.m_hwm, nrows);
    counts = nrows;
    for (auto &v: counts) v *= m_row_dsize;
    mpi::counts_to_displs_consec(counts, displs);
    auto nrow_total = std::accumulate(nrows.cbegin(), nrows.cend(), 0ul);
    push_back(nrow_total);
    mpi::all_gatherv(src.dbegin(), m_hwm * m_row_dsize, dbegin(), counts, displs);
    post_insert_range(0, nrow_total);
}

void TableBase::erase_protector(RowProtector *rp) const {
    m_row_protectors.erase(rp->m_it);
}

bool TableBase::is_protected() const {
    size_t protected_here = std::any_of(m_row_protectors.cbegin(), m_row_protectors.cend(),
                       [](const RowProtector* rp){return rp->is_protected();});
    return mpi::all_lor(protected_here);
}

bool TableBase::is_protected(const size_t& irow) const {
    return std::any_of(m_row_protectors.cbegin(), m_row_protectors.cend(),
                       [&](const RowProtector* rp){return rp->is_protected(irow);});
}

TableBase::Loc::Loc(size_t irank, size_t irow) : m_irank(irank), m_irow(irow){
#ifndef DNDEBUG
    mpi::bcast(irank);
    mpi::bcast(irow);
    MPI_ASSERT_ALL(m_irank==irank, "rank index in TableBase::Loc should be consistent across all ranks");
    MPI_ASSERT_ALL(m_irow==irow, "row index in TableBase::Loc should be consistent across all ranks");
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
