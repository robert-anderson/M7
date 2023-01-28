//
// Created by Robert J. Anderson on 08/05/2021.
//

#include <numeric>
#include <M7_lib/io/Logging.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/sort/ExtremalIndices.h>

#include "TableBase.h"
#include "RowProtector.h"
#include "M7_lib/util/Pointer.h"

void TableBase::clear(uint_t i) {
    DEBUG_ASSERT_LT(i, nrow_in_use(), "row index OOB");
    DEBUG_ASSERT_FALSE(is_protected(i), "cannot clear a protected row");
    if (m_bw.i_can_modify()) std::memset(begin(i), 0, row_size());
}

bool TableBase::is_clear(uint_t i) const {
    return std::memcmp(cbegin(i), m_null_row_string.data(), row_size()) == 0;
}

TableBase::TableBase(uint_t row_size) : m_bw(row_size), m_null_row_string(row_size, 0){}

TableBase::TableBase(const TableBase &other) : TableBase(other.m_bw.m_row_size){}

void TableBase::clear() {
    REQUIRE_FALSE(is_protected(), "cannot clear a table with protected records");
    m_bw.clear();
}

bool TableBase::is_clear() const {
    for (uint_t irow=0ul; irow<m_bw.m_nrow; ++irow) if (!is_clear(irow)) return false;
    return true;
}

void TableBase::set_buffer(Buffer *buffer) {
    DEBUG_ASSERT_TRUE(buffer, "buffer is null");
    DEBUG_ASSERT_TRUE(!m_bw.allocated(), "buffer window is not allocated")
    buffer->append_window(&m_bw);
}

uint_t TableBase::push_back(uint_t n) {
    DEBUG_ASSERT_TRUE(n, "cannot push_back by zero rows");
    DEBUG_ASSERT_TRUE(row_size(), "cannot resize a table with zero record size");
    const auto nbyte_add = n * m_bw.m_row_size;
    if (!m_bw.cend() || ptr::after_end(m_bw.cend() + nbyte_add, m_bw.cbegin()+m_bw.m_size)) expand(n);
    const auto tmp = nrow_in_use();
    m_bw.set_end(m_bw.m_nrow + n);
    return tmp;
}

uint_t TableBase::get_free_row() {
    if (m_freed_rows.empty()) return push_back();
    const auto irec = m_freed_rows.top();
    m_freed_rows.pop();
    m_is_freed_row[irec] = false;
    return irec;
}

void TableBase::free(uint_t i) {
    DEBUG_ASSERT_LT(i, nrow_in_use(), "record index OOB");
    DEBUG_ASSERT_FALSE(is_protected(i), "cannot clear a protected record");
    DEBUG_ASSERT_FALSE(is_freed(i), "cannot free a record that is already freed");
    clear(i);
    m_freed_rows.push(i);
    m_is_freed_row[i] = true;
}

void TableBase::free() {
    clear();
    while (!m_freed_rows.empty()) m_freed_rows.pop();
    m_is_freed_row.assign(m_is_freed_row.size(), false);
}

bool TableBase::empty() const {
    return m_bw.empty();
}

bool TableBase::is_freed(uint_t i) const {
    DEBUG_ASSERT_LT(i, nrow_in_use(), "record index OOB");
    return m_is_freed_row[i];
}

uint_t TableBase::bw_size() const {
    return m_bw.m_size;
}

void TableBase::resize(uint_t nrec, double factor) {
    REQUIRE_TRUE(nrec, "new size should be non-zero");
    REQUIRE_TRUE(row_size(), "cannot resize, row size is zero");
    REQUIRE_GE(nrec, nrow_in_use(), "resize would discard uncleared data");
    m_bw.resize(nrec * row_size(), factor);
    DEBUG_ASSERT_LT(this->nrow_in_use(), m_bw.m_size / row_size(), "resize has discarded uncleared data");
    m_is_freed_row.resize(capacity(), false);
}

void TableBase::expand(uint_t nrow, double factor) {
    resize(this->capacity() + nrow, factor);
}

void TableBase::free_many(const uintv_t &irows) {
    for (auto irec : irows) free(irec);
}

#if 0
void TableBase::insert_records(const Buffer::Window &recv, uint_t nrec, const std::list<recv_cb_t> &callbacks) {
    for (uint_t irec_recv = 0; irec_recv < nrec; ++irec_recv) {
        auto irec = get_free_row();
        std::memcpy(begin(irec), recv.m_begin + irec_recv * row_size(), row_size());
        post_insert(irec);
        for (auto f: callbacks) f(irec);
    }
}
void TableBase::transfer_records(const uintv_t &/*irecs*/, uint_t /*irank_send*/,
                                 uint_t /*irank_recv*/, const std::list<recv_cb_t>& /*callbacks*/){
    DEBUG_ASSERT_NE_ALL(irank_recv, irank_send, "sending and recving ranks should never be the same");
    if (!m_transfer) m_transfer = ptr::smart::make_unique<RowTransfer>(m_bw.name());
    uint_t nrec = 0;
    if (mpi::i_am(irank_send)){
        auto& send_bw = m_transfer->m_send_bw;
        nrec = irecs.size();
        mpi::send(&nrec, 1, irank_recv, m_transfer->m_nrec_p2p_tag);
        if (!nrec){
            logging::debug_("Sending rank notifying recving rank that no records are transferred");
            return;
        }
        logging::info_("Transferring {} records outward to rank {}", nrec, irank_recv);

        auto size_required = nrec * row_size();
        if (send_bw.m_size < size_required) send_bw.resize(size_required);
        for (auto iirec = 0ul; iirec < nrec; ++iirec) {
            const auto &irec = irecs[iirec];
            std::memcpy(send_bw.m_begin + iirec * row_size(), begin(irec), row_size());
        }
        mpi::send(send_bw.m_begin, row_size() * nrec, irank_recv, m_transfer->m_irecs_p2p_tag);
        /*
         * sent records can now be erased
         */
        free_many(irecs);
    }
    if (mpi::i_am(irank_recv)){
        auto& recv_bw = m_transfer->m_send_bw;
        mpi::recv(&nrec, 1, irank_send, m_transfer->m_nrec_p2p_tag);
        if (!nrec){
            logging::debug_("Recving rank notified by sending rank that no records are transferred");
            return;
        }
        auto size_required = nrec * row_size();
        if (recv_bw.m_size<size_required) recv_bw.resize(size_required);
        logging::info_("Transferring {} records inward from rank {}", nrec, irank_send);
        mpi::recv(recv_bw.m_begin, row_size() * nrec, irank_send, m_transfer->m_irecs_p2p_tag);
        /*
         * now emplace received records in TableBase buffer window, and call all callbacks for each
         */
        insert_records(recv_bw, nrec, callbacks);
    }
}
#endif

void TableBase::copy_record_in(const TableBase &src, uint_t isrc, uint_t idst) {
    DEBUG_ASSERT_LT(isrc, src.nrow_in_use(), "src record index OOB");
    DEBUG_ASSERT_LT(idst, nrow_in_use(), "dst record index OOB");
    if (!m_bw.i_can_modify()) return;
    std::memcpy(begin(idst), src.cbegin(isrc), row_size());
}

void TableBase::swap_records(uint_t i, uint_t j) {
    if (!m_bw.i_can_modify()) return;
    if (i == j) return;
    auto iptr = begin(i);
    auto jptr = begin(j);
    std::swap_ranges(iptr, iptr + row_size(), jptr);
}

str_t TableBase::to_string(const uintv_t *ordering) const {
    str_t out;
    auto begin_ptr = cbegin();
    for (uint_t i=0ul; i<nrow_in_use(); ++i){
        auto irec = ordering ? ordering->at(i) : i;
        auto rec = begin_ptr + irec * row_size();
        for (uint_t ibyte=0ul; ibyte < row_size(); ++ibyte){
            out+= std::to_string(static_cast<int>(rec[ibyte])) + " ";
        }
        out+="\n";
    }
    return out;
}

void TableBase::all_gatherv(const TableBase &src) {
    TableBase::clear();
    uintv_t nrecs(mpi::nrank());
    uintv_t counts(mpi::nrank());
    uintv_t displs(mpi::nrank());
    REQUIRE_EQ_ALL(src.row_size(), row_size(),
                    "the size of records being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.nrow_in_use(), nrecs);
    counts = nrecs;
    for (auto &v: counts) v *= row_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrec_total = std::accumulate(nrecs.cbegin(), nrecs.cend(), 0ul);
    if (!nrec_total) return;
    push_back(nrec_total);
    if (m_bw.node_shared()) {
        /*
         * if a subset of rows (just the node roots) are gathering, we can't use gatherv (one recving rank) or
         * all_gatherv (all ranks gather), so we have to use the general all_to_allv
         */
        // initially, assume all ranks receive all data from this rank
        uintv_t sendcounts(mpi::nrank(), src.nrow_in_use()*row_size());
        // set to zero all sendcounts for ranks that are not node-roots
        for (uint_t irank=0ul; irank<mpi::nrank(); ++irank){
            if (!mpi::is_node_root(irank)) sendcounts[irank] = 0;
        }
        // all send displacements remain at 0ul
        uintv_t senddispls(mpi::nrank(), 0ul);
        uintv_t recvcounts(mpi::nrank(), 0ul);
        uintv_t recvdispls(mpi::nrank(), 0ul);
        // only node-roots receive non-zero counts from senders
        if (mpi::on_node_i_am_root()) recvcounts = counts;
        recvdispls = mpi::counts_to_displs_consec(recvcounts);
        mpi::all_to_allv(src.cbegin(), sendcounts, senddispls, begin(), recvcounts, recvdispls);
    }
    else {
        mpi::all_gatherv(cbegin(), src.m_bw.size_in_use(), begin(), counts, displs);
    }
}

void TableBase::gatherv(const TableBase &src, uint_t irank) {
    if (mpi::i_am(irank)) {
        REQUIRE_TRUE(m_bw.i_can_modify(),
                     "node-shared memory: cannot gather to a rank which is not the root on its node");
        TableBase::clear();
    }
    uintv_t nrecs(mpi::nrank());
    uintv_t counts(mpi::nrank());
    uintv_t displs(mpi::nrank());
    REQUIRE_EQ_ALL(src.row_size(), row_size(),
                    "the size of records being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.nrow_in_use(), nrecs);
    counts = nrecs;
    for (auto &v: counts) v *= row_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrec_total = std::accumulate(nrecs.cbegin(), nrecs.cend(), 0ul);
    if (mpi::i_am(irank)) push_back(nrec_total);
    mpi::gatherv(src.cbegin(), src.m_bw.size_in_use(), begin(), counts, displs, irank);
}

bool TableBase::is_protected() const {
    return !m_protected_rows.empty();
}

bool TableBase::freed_rows_consistent() const {
    const auto nfree_vec = std::accumulate(m_is_freed_row.cbegin(), m_is_freed_row.cend(), 0ul);
    if (nfree_vec != m_freed_rows.size()) return false;
    auto tmp_stack = m_freed_rows;
    while (!tmp_stack.empty()) {
        if (!m_is_freed_row[tmp_stack.top()]) return false;
        tmp_stack.pop();
    }
    return true;
}

TableBase::Loc::Loc(uint_t irank, uint_t irec) : m_irank(irank), m_irec(irec){}

TableBase::Loc::Loc(uint_t irec) : Loc(mpi::all_min(uint_t(irec==~0ul ? mpi::nrank() : mpi::irank())), irec) {
    const uint_t is_mine = (irec != ~0ul);
    REQUIRE_EQ_ALL(mpi::all_sum(is_mine), 1ul,
                   "a distributed table location can only be on one rank!");
}

TableBase::Loc::operator bool() const {
    return m_irank!=~0ul;
}

bool TableBase::Loc::is_mine() const {
    return mpi::i_am(m_irank);
}

bool TableBase::Loc::operator==(const TableBase::Loc &other) {
    return m_irank==other.m_irank and m_irec == other.m_irec;
}

bool TableBase::Loc::operator!=(const TableBase::Loc &other) {
    return !(*this==other);
}

