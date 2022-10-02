//
// Created by Robert J. Anderson on 08/05/2021.
//

#include <numeric>
#include <M7_lib/io/Logging.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/sort/ExtremalIndices.h>

#include "TableBase.h"
#include "RowProtector.h"
#include "M7_lib/util/SmartPtr.h"

TableBase::TableBase(uint_t record_size) :
        m_bw(record_size), m_null_record_string(record_size, 0){}

TableBase::TableBase(const TableBase &other) : TableBase(other.m_bw.m_record_size){}

buf_t *TableBase::begin() {
    return m_bw.m_begin;
}

const buf_t *TableBase::begin() const {
    return m_bw.m_begin;
}

buf_t *TableBase::begin(uint_t irec) {
    return m_bw.m_begin + irec * record_size();
}

const buf_t *TableBase::begin(uint_t irec) const {
    return m_bw.m_begin + irec * record_size();
}

void TableBase::set_buffer(Buffer *buffer) {
    ASSERT(buffer);
    ASSERT(!m_bw.allocated())
    buffer->append_window(&m_bw);
}

bool TableBase::is_full() const {
    return m_hwm == nrecord();
}

uint_t TableBase::push_back(uint_t nrec) {
    DEBUG_ASSERT_TRUE(record_size(), "cannot resize a table with zero record size");
    if (m_hwm + nrec > this->nrecord()) expand(nrec);
    auto tmp = m_hwm;
    m_hwm += nrec;
    return tmp;
}

uint_t TableBase::get_empty_record() {
    if (m_empty_records.empty()) return push_back();
    auto irec = m_empty_records.top();
    m_empty_records.pop();
    return irec;
}

void TableBase::clear() {
    ASSERT(!is_protected());
    if (!m_bw.allocated()) return;
    std::memset(begin(), 0, record_size() * m_hwm);
    m_hwm = 0ul;
    while (!m_empty_records.empty()) m_empty_records.pop();
}

void TableBase::clear(uint_t irec) {
    DEBUG_ASSERT_LT(irec, m_hwm, "row index OOB");
    DEBUG_ASSERT_FALSE(is_protected(irec), "cannot clear a protected row");
    std::memset(begin(irec), 0, record_size());
    m_empty_records.push(irec);
}

bool TableBase::is_cleared() const {
    return std::memcmp(begin(), m_null_record_string.data(), record_size()) == 0;
}

bool TableBase::is_cleared(uint_t irec) const {
    return std::memcmp(begin(irec), m_null_record_string.data(), record_size()) == 0;
}

uint_t TableBase::bw_size() const {
    return m_bw.m_size;
}

void TableBase::resize(uint_t nrec, double factor) {
    DEBUG_ASSERT_TRUE(record_size(), "cannot resize, row size is zero");
    DEBUG_ASSERT_GE(nrec, m_hwm, "resize would discard uncleared data");
    m_bw.resize(nrec * record_size(), factor);
    DEBUG_ASSERT_LT(m_hwm, m_bw.m_size / record_size(), "resize has discarded uncleared data");
}

void TableBase::expand(uint_t nrec, double factor) {
    resize(this->nrecord() + nrec, factor);
}

void TableBase::clear_records(const uintv_t &irecs) {
    for (auto irec : irecs) clear(irec);
}

void TableBase::insert_records(const Buffer::Window &recv, uint_t nrec, const std::list<recv_cb_t> &callbacks) {
    for (uint_t irec_recv = 0; irec_recv < nrec; ++irec_recv) {
        auto irec = get_empty_record();
        std::memcpy(begin(irec), recv.m_begin + irec_recv * record_size(), record_size());
        post_insert(irec);
        for (auto f: callbacks) f(irec);
    }
}
void TableBase::transfer_records(const uintv_t &irecs, uint_t irank_send,
                                 uint_t irank_recv, const std::list<recv_cb_t>& callbacks){
    DEBUG_ASSERT_NE_ALL(irank_recv, irank_send, "sending and recving ranks should never be the same");
    if (!m_transfer) m_transfer = smart_ptr::make_unique<RowTransfer>(m_bw.name());
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

        auto size_required = nrec * record_size();
        if (send_bw.m_size < size_required) send_bw.resize(size_required);
        for (auto iirec = 0ul; iirec < nrec; ++iirec) {
            const auto &irec = irecs[iirec];
            std::memcpy(send_bw.m_begin + iirec * record_size(), begin(irec), record_size());
        }
        mpi::send(send_bw.m_begin, record_size() * nrec, irank_recv, m_transfer->m_irecs_p2p_tag);
        /*
         * sent records can now be erased
         */
        clear_records(irecs);
    }
    if (mpi::i_am(irank_recv)){
        auto& recv_bw = m_transfer->m_send_bw;
        mpi::recv(&nrec, 1, irank_send, m_transfer->m_nrec_p2p_tag);
        if (!nrec){
            logging::debug_("Recving rank notified by sending rank that no records are transferred");
            return;
        }
        auto size_required = nrec * record_size();
        if (recv_bw.m_size<size_required) recv_bw.resize(size_required);
        logging::info_("Transferring {} records inward from rank {}", nrec, irank_send);
        mpi::recv(recv_bw.m_begin, record_size() * nrec, irank_send, m_transfer->m_irecs_p2p_tag);
        /*
         * now emplace received records in TableBase buffer window, and call all callbacks for each
         */
        insert_records(recv_bw, nrec, callbacks);
    }
}

void TableBase::copy_record_in(const TableBase &src, uint_t irec_src, uint_t irec_dst) {
    ASSERT(irec_dst < m_hwm);
    std::memcpy(begin(irec_dst), src.begin(irec_src), record_size());
}

void TableBase::swap_records(uint_t irec, uint_t jrec) {
    if (irec == jrec) return;
    auto iptr = begin(irec);
    auto jptr = begin(jrec);
    std::swap_ranges(iptr, iptr + record_size(), jptr);
}

str_t TableBase::to_string(const uintv_t *ordering) const {
    str_t out;
    auto begin_ptr = begin();
    for (uint_t i=0ul; i<m_hwm; ++i){
        auto irec = ordering ? ordering->at(i) : i;
        auto rec = begin_ptr + irec * record_size();
        for (uint_t ibyte=0ul; ibyte < record_size(); ++ibyte){
            out+= std::to_string(static_cast<int>(rec[ibyte])) + " ";
        }
        out+="\n";
    }
    return out;
}

void TableBase::all_gatherv(const TableBase &src) {
    clear();
    uintv_t nrecs(mpi::nrank());
    uintv_t counts(mpi::nrank());
    uintv_t displs(mpi::nrank());
    DEBUG_ASSERT_EQ(src.record_size(), record_size(),
                    "the size of records being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.m_hwm, nrecs);
    counts = nrecs;
    for (auto &v: counts) v *= record_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrec_total = std::accumulate(nrecs.cbegin(), nrecs.cend(), 0ul);
    push_back(nrec_total);
    mpi::all_gatherv(src.begin(), src.m_hwm * record_size(), begin(), counts, displs);
    post_insert_range(0, nrec_total);
}

void TableBase::gatherv(const TableBase &src, uint_t irank) {
    if (mpi::i_am(irank)) clear();
    uintv_t nrecs(mpi::nrank());
    uintv_t counts(mpi::nrank());
    uintv_t displs(mpi::nrank());
    DEBUG_ASSERT_EQ(src.record_size(), record_size(),
                    "the size of records being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.m_hwm, nrecs);
    counts = nrecs;
    for (auto &v: counts) v *= record_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrec_total = std::accumulate(nrecs.cbegin(), nrecs.cend(), 0ul);
    if (mpi::i_am(irank))
        push_back(nrec_total);
    mpi::gatherv(src.begin(), src.m_hwm * record_size(), begin(), counts, displs, irank);
    if (mpi::i_am(irank))
        post_insert_range(0, nrec_total);
}

bool TableBase::is_protected() const {
    return !m_protected_records.empty();
}

uint_t TableBase::nrecord_nonempty() const {
    return m_hwm - m_empty_records.size();
}

std::set<uint_t> TableBase::empty_records_set() const {
    std::set<uint_t> set;
    auto stack = m_empty_records;
    while (!stack.empty()) {
        set.insert(stack.top());
        stack.pop();
    }
    return set;
}

TableBase::Loc::Loc(uint_t irank, uint_t irec) : m_irank(irank), m_irec(irec){
#ifndef NDEBUG
    mpi::bcast(irank);
    mpi::bcast(irec);
    DEBUG_ASSERT_EQ(m_irank, irank, "rank index in TableBase::Loc should be consistent across all ranks");
    DEBUG_ASSERT_EQ(m_irec, irec, "record index in TableBase::Loc should be consistent across all ranks");
#endif
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
