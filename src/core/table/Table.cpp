//
// Created by rja on 09/02/2021.
//

#include "Table.h"
#include "src/core/io/Logging.h"
#include "src/core/parallel/MPIAssert.h"
#include "src/core/sort/ExtremalValues.h"


TableBase::TableBase(size_t row_dsize) :
        m_row_dsize(row_dsize), m_row_size(row_dsize * defs::nbyte_data){}

TableBase::TableBase(const TableBase &other) :
        TableBase(other.m_row_dsize){}

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
    if (!m_bw.allocated()) return;
    std::memset(dbegin(), 0, m_row_size * m_hwm);
    m_hwm = 0ul;
    while (!m_free_rows.empty()) m_free_rows.pop();
}

void TableBase::clear(const size_t &irow) {
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

/*
void TableBase::print_contents(const defs::inds *ordering) const {
    const auto n = ordering ? std::min(ordering->size(), m_hwm) : m_hwm;
    for (size_t iirow = 0ul; iirow < n; ++iirow) {
        auto irow = ordering ? (*ordering)[iirow] : iirow;
        std::cout << irow << ". ";
        for (auto column: m_columns) {
            std::cout << column->to_string(irow) + " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void TableBase::print_contents(const ExtremalValues &xv) const {
    defs::inds tmp;
    tmp.reserve(xv.nfound());
    for (size_t i = 0ul; i < xv.nfound(); ++i) tmp.push_back(xv[i]);
    print_contents(&tmp);
}
*/

void TableBase::resize(size_t nrow) {
    assert(nrow > m_nrow);
    m_bw.resize(nrow * m_row_dsize);
    m_nrow = nrow;
}

void TableBase::expand(size_t nrow, double expansion_factor) {
    m_bw.expand(nrow * m_row_dsize, expansion_factor);
    m_nrow = m_bw.dsize()/m_row_dsize;
}

void TableBase::expand(size_t nrow) {
    m_bw.expand(nrow * m_row_dsize);
    m_nrow = m_bw.dsize()/m_row_dsize;
}

void TableBase::erase_rows(const defs::inds &irows) {
    for (auto irow : irows) {
        clear(irow);
    }
}

void TableBase::post_insert(const size_t& iinsert) {}

void TableBase::insert_rows(const Buffer::Window &recv, size_t nrow, const std::list<recv_cb_t> &callbacks) {
    for (size_t irow_recv = 0; irow_recv < nrow; ++irow_recv) {
        auto irow_TableX = get_free_row();
        std::memcpy(dbegin(irow_TableX), recv.m_dbegin + irow_recv * m_row_dsize, m_row_size);
        post_insert(irow_TableX);
        for (auto f: callbacks) f(irow_TableX);
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
