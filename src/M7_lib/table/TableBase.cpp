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

TableBase::TableBase(uint_t slot_size) : m_bw(slot_size){}

TableBase::TableBase(const TableBase &other) : TableBase(other.m_bw.m_slot_size){}

buf_t *TableBase::begin() {
    return m_bw.m_begin;
}

const buf_t *TableBase::begin() const {
    return m_bw.m_begin;
}

buf_t *TableBase::begin(uint_t islot) {
    return m_bw.m_begin + islot * slot_size();
}

const buf_t *TableBase::begin(uint_t islot) const {
    return m_bw.m_begin + islot * slot_size();
}

void TableBase::set_buffer(Buffer *buffer) {
    ASSERT(buffer);
    ASSERT(!m_bw.allocated())
    buffer->append_window(&m_bw);
}

uint_t TableBase::push_back(uint_t nslot) {
    DEBUG_ASSERT_TRUE(slot_size(), "cannot resize a table with zero record size");
    if (m_hwm + nslot > this->nslot()) expand(nslot);
    auto tmp = m_hwm;
    m_hwm += nslot;
    return tmp;
}

uint_t TableBase::get_free_slot() {
    if (m_freed_slots.empty()) return push_back();
    auto irec = m_freed_slots.top();
    m_freed_slots.pop();
    return irec;
}

void TableBase::clear() {
    DEBUG_ASSERT_FALSE(is_protected(), "cannot clear a table with protected records");
    if (!m_bw.allocated()) return;
    std::memset(begin(), 0, slot_size() * m_hwm);
    m_hwm = 0ul;
    while (!m_freed_slots.empty()) m_freed_slots.pop();
    m_is_freed_slot.assign(m_is_freed_slot.size(), false);
}

void TableBase::free(uint_t islot) {
    DEBUG_ASSERT_LT(islot, m_hwm, "slot index OOB");
    DEBUG_ASSERT_FALSE(is_protected(islot), "cannot clear a protected record");
    DEBUG_ASSERT_FALSE(is_freed(islot), "cannot free a slot that is already freed");
    std::memset(begin(islot), 0, slot_size());
    m_freed_slots.push(islot);
    m_is_freed_slot[islot] = true;
}

bool TableBase::empty() const {
    return !m_hwm;
}

bool TableBase::is_freed(uint_t islot) const {
    DEBUG_ASSERT_LT(islot, m_hwm, "slot index OOB");
    return m_is_freed_slot[islot];
}

uint_t TableBase::bw_size() const {
    return m_bw.m_size;
}

void TableBase::resize(uint_t nrec, double factor) {
    DEBUG_ASSERT_TRUE(nrec, "new size should be non-zero");
    DEBUG_ASSERT_TRUE(slot_size(), "cannot resize, row size is zero");
    DEBUG_ASSERT_GE(nrec, m_hwm, "resize would discard uncleared data");
    m_bw.resize(nrec * slot_size(), factor);
    DEBUG_ASSERT_LT(m_hwm, m_bw.m_size / slot_size(), "resize has discarded uncleared data");
    m_is_freed_slot.resize(nslot(), false);
}

void TableBase::expand(uint_t nrec, double factor) {
    resize(this->nslot() + nrec, factor);
}

void TableBase::free_records(const uintv_t &irecs) {
    for (auto irec : irecs) free(irec);
}

void TableBase::insert_records(const Buffer::Window &recv, uint_t nrec, const std::list<recv_cb_t> &callbacks) {
    for (uint_t irec_recv = 0; irec_recv < nrec; ++irec_recv) {
        auto irec = get_free_slot();
        std::memcpy(begin(irec), recv.m_begin + irec_recv * slot_size(), slot_size());
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

        auto size_required = nrec * slot_size();
        if (send_bw.m_size < size_required) send_bw.resize(size_required);
        for (auto iirec = 0ul; iirec < nrec; ++iirec) {
            const auto &irec = irecs[iirec];
            std::memcpy(send_bw.m_begin + iirec * slot_size(), begin(irec), slot_size());
        }
        mpi::send(send_bw.m_begin, slot_size() * nrec, irank_recv, m_transfer->m_irecs_p2p_tag);
        /*
         * sent records can now be erased
         */
        free_records(irecs);
    }
    if (mpi::i_am(irank_recv)){
        auto& recv_bw = m_transfer->m_send_bw;
        mpi::recv(&nrec, 1, irank_send, m_transfer->m_nrec_p2p_tag);
        if (!nrec){
            logging::debug_("Recving rank notified by sending rank that no records are transferred");
            return;
        }
        auto size_required = nrec * slot_size();
        if (recv_bw.m_size<size_required) recv_bw.resize(size_required);
        logging::info_("Transferring {} records inward from rank {}", nrec, irank_send);
        mpi::recv(recv_bw.m_begin, slot_size() * nrec, irank_send, m_transfer->m_irecs_p2p_tag);
        /*
         * now emplace received records in TableBase buffer window, and call all callbacks for each
         */
        insert_records(recv_bw, nrec, callbacks);
    }
}

void TableBase::copy_record_in(const TableBase &src, uint_t islot_src, uint_t islot_dst) {
    ASSERT(islot_dst < m_hwm);
    std::memcpy(begin(islot_dst), src.begin(islot_src), slot_size());
}

void TableBase::swap_records(uint_t irec, uint_t jrec) {
    if (irec == jrec) return;
    auto iptr = begin(irec);
    auto jptr = begin(jrec);
    std::swap_ranges(iptr, iptr + slot_size(), jptr);
}

str_t TableBase::to_string(const uintv_t *ordering) const {
    str_t out;
    auto begin_ptr = begin();
    for (uint_t i=0ul; i<m_hwm; ++i){
        auto irec = ordering ? ordering->at(i) : i;
        auto rec = begin_ptr + irec * slot_size();
        for (uint_t ibyte=0ul; ibyte < slot_size(); ++ibyte){
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
    DEBUG_ASSERT_EQ(src.slot_size(), slot_size(),
                    "the size of records being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.m_hwm, nrecs);
    counts = nrecs;
    for (auto &v: counts) v *= slot_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrec_total = std::accumulate(nrecs.cbegin(), nrecs.cend(), 0ul);
    push_back(nrec_total);
    mpi::all_gatherv(src.begin(), src.m_hwm * slot_size(), begin(), counts, displs);
    post_insert_range(0, nrec_total);
}

void TableBase::gatherv(const TableBase &src, uint_t irank) {
    if (mpi::i_am(irank)) clear();
    uintv_t nrecs(mpi::nrank());
    uintv_t counts(mpi::nrank());
    uintv_t displs(mpi::nrank());
    DEBUG_ASSERT_EQ(src.slot_size(), slot_size(),
                    "the size of records being gathered does not match that stored in the gathering table");
    mpi::all_gather(src.m_hwm, nrecs);
    counts = nrecs;
    for (auto &v: counts) v *= slot_size();
    mpi::counts_to_displs_consec(counts, displs);
    auto nrec_total = std::accumulate(nrecs.cbegin(), nrecs.cend(), 0ul);
    if (mpi::i_am(irank))
        push_back(nrec_total);
    mpi::gatherv(src.begin(), src.m_hwm * slot_size(), begin(), counts, displs, irank);
    if (mpi::i_am(irank))
        post_insert_range(0, nrec_total);
}

bool TableBase::is_protected() const {
    return !m_protected_records.empty();
}

uint_t TableBase::nrecord() const {
    return m_hwm - m_freed_slots.size();
}

bool TableBase::freed_slots_consistent() const {
    const auto nfree_vec = std::accumulate(m_is_freed_slot.cbegin(), m_is_freed_slot.cend(), 0ul);
    if (nfree_vec != m_freed_slots.size()) return false;
    auto tmp_stack = m_freed_slots;
    while (!tmp_stack.empty()) {
        if (!m_is_freed_slot[tmp_stack.top()]) return false;
        tmp_stack.pop();
    }
    return true;
}

TableBase::Loc::Loc(uint_t irank, uint_t islot) : m_irank(irank), m_islot(islot){
#ifndef NDEBUG
    mpi::bcast(irank);
    mpi::bcast(islot);
    DEBUG_ASSERT_EQ(m_irank, irank, "rank index in TableBase::Loc should be consistent across all ranks");
    DEBUG_ASSERT_EQ(m_islot, islot, "record index in TableBase::Loc should be consistent across all ranks");
#endif
}

TableBase::Loc::operator bool() const {
    return m_irank!=~0ul;
}

bool TableBase::Loc::is_mine() const {
    return mpi::i_am(m_irank);
}

bool TableBase::Loc::operator==(const TableBase::Loc &other) {
    return m_irank==other.m_irank and m_islot == other.m_islot;
}

bool TableBase::Loc::operator!=(const TableBase::Loc &other) {
    return !(*this==other);
}
