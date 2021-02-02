//
// Created by rja on 02/02/2021.
//

#include "RankAllocator.h"

void RankAllocatorBase::refresh_callback_list() {
    m_recv_callbacks.clear();
    for (auto ptr : m_dependents) {
        m_recv_callbacks.emplace_back([ptr](size_t irow) { ptr->on_row_recv_(irow); });
    }
}

std::list<RankAllocatorBase::Dynamic *>::iterator
RankAllocatorBase::add_dependent(RankAllocatorBase::Dynamic *dependent) {
    m_dependents.push_back(dependent);
    auto it = m_dependents.end();
    refresh_callback_list();
    return --it;
}

void RankAllocatorBase::erase_dependent(RankAllocatorBase::Dynamic *dependent) {
    m_dependents.erase(dependent->m_it);
    refresh_callback_list();
}

RankAllocatorBase::RankAllocatorBase(Table &table, size_t nblock, size_t period) :
        m_table(table), m_nblock(nblock), m_period(period),
        m_block_to_rank(nblock, 0ul), m_rank_to_blocks(mpi::nrank()),
        m_mean_work_times(nblock, 0.0), m_gathered_times(mpi::nrank(), 0.0)
{
    MPI_REQUIRE_ALL(m_acceptable_imbalance>=0.0, "Acceptable imbalance fraction must be non-negative");
    MPI_REQUIRE_ALL(m_acceptable_imbalance<=1.0, "Acceptable imbalance fraction must not exceed 100%");
    size_t iblock = 0ul;
    for (auto &rank:m_block_to_rank) {
        rank = iblock%mpi::nrank();
        m_rank_to_blocks[rank].push_front(iblock);
        iblock++;
    }
}

size_t RankAllocatorBase::nblock_() const {
    return m_rank_to_blocks[mpi::irank()].size();
}

void RankAllocatorBase::record_work_time(const size_t &irow, const Timer &work_time) {
    m_mean_work_times[get_block_irow(irow)]+=work_time;
}

size_t RankAllocatorBase::get_nskip_() const {
    /*
     * only called on the sending rank since the block-wise work times are not
     * shared with all ranks.
     *
     * When selecting a block to transfer, we don't want to send an expensive block
     * as this has the potential to be counterproductive.
     *
     * Thus, send the first block found by iterating from the front of the list with
     * mean_work_time less than the average over all blocks on this MPI rank.
     */

    const auto& times = m_mean_work_times;
    const auto& list = m_rank_to_blocks[mpi::irank()];

    double mean = 0.0;
    for (auto it=list.cbegin(); it!=list.cend(); ++it) mean+=times[*it];
    mean/=list.size();
    size_t nskip=0ul;
    for (auto it=list.cbegin(); it!=list.cend(); ++it){
        if (times[*it] < mean) return nskip;
        ++nskip;
    }
    /*
     * if the method has not already returned, all blocks have the same mean_work_time.
     */
    return 0;
}

void RankAllocatorBase::deactivate() {
    log::info("Deactivating dynamic load balancing.");
    m_icycle_active = ~0ul;
}

void RankAllocatorBase::activate(size_t icycle) {
    log::info("Activating dynamic load balancing.");
    m_icycle_active = icycle;
    // rezero the counter in case of later reactivation
    m_nnull_updates = 0ul;
    m_mean_work_times.assign(m_nblock, 0.0);
}

bool RankAllocatorBase::is_active() const {
    return m_icycle_active!=~0ul;
}

bool RankAllocatorBase::consistent() {
    for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
        for (auto& iblock: m_rank_to_blocks[irank]){
            if (m_block_to_rank[iblock]!=irank) return false;
        }
    }
    return true;
}

void RankAllocatorBase::update(size_t icycle) {
    if (!is_active()) return;
    size_t ncycle_active = m_icycle_active-icycle;
    if (!ncycle_active || ncycle_active%m_period) return;
    // else, this is a preparatory iteration
    auto local_time = std::accumulate(m_mean_work_times.begin(), m_mean_work_times.end(), 0.0);
    mpi::all_gather(local_time, m_gathered_times);

    auto it_lazy = std::min_element(m_gathered_times.begin(), m_gathered_times.end());
    auto it_busy = std::max_element(m_gathered_times.begin(), m_gathered_times.end());
    if (*it_lazy==0.0) {
        log::warn("The most idle rank appears to have done no work at all");
        log::debug("Gathered times: {}", utils::to_string(m_gathered_times));
    }
    MPI_REQUIRE_ALL(*it_busy>0.0, "The busiest rank appears to have done no work. "
                                  "record_work_time must be called by application");

    // this rank has done the least work, so it should receive a block
    size_t irank_recv = std::distance(m_gathered_times.begin(), it_lazy);
    // this rank has done the most work, so it should send a block
    size_t irank_send = std::distance(m_gathered_times.begin(), it_busy);
    if (m_rank_to_blocks[irank_send].size()==1){
        log::warn("Busiest rank has only one (very expensive) block remaining ({})",
                  m_rank_to_blocks[irank_send].front());
        log::warn("Load balance cannot be further improved");
        /*
         * if the cause is many expensive rows in a single block, they can probably be broken up by
         * increasing the number of blocks. otherwise, there's nothing to be done to improve parallel
         * performance since the system in question cannot be balanced its current representation.
         */
        deactivate();
        return;
    }

    if (irank_send==irank_recv) return;
    auto realloc_ratio = 1-m_acceptable_imbalance;
    if (m_gathered_times[irank_recv] > realloc_ratio * m_gathered_times[irank_send]) {
        ++m_nnull_updates;
        if(m_nnull_updates>=m_nnull_updates_deactivate) {
            log::info(
                    "Load imbalance has been below {}% for over {} periods of {} cycles.",
                    m_acceptable_imbalance * 100, m_nnull_updates, m_period);
            deactivate();
        }
        m_mean_work_times.assign(m_nblock, 0.0);
        return;
    }
    else {
        // not a null update: rezero the counter
        m_nnull_updates = 0ul;
    }
    MPI_REQUIRE_ALL(!m_rank_to_blocks[irank_send].empty(),
                    "Sending rank has no blocks to send! Specified work_time value may be erroneous");
    /*
     * sender will send all rows in the selected block, which will become the front
     * block of the recver.
     */
    size_t nskip;
    if (mpi::i_am(irank_send)) nskip = get_nskip_();
    mpi::bcast(nskip, irank_send);
    MPI_ASSERT(nskip<m_rank_to_blocks[irank_send].size(), "Too many skips");
    auto it_block_transfer = m_rank_to_blocks[irank_send].begin();
    for (size_t iskip=0ul; iskip<nskip; ++iskip) ++it_block_transfer;

    log::info("Sending block {} from rank {} to {} on cycle {}", *it_block_transfer, irank_send, irank_recv, icycle);

    // prepare vector of row indices to send
    defs::inds irows_send;
    if (mpi::i_am(irank_send)){
        for (size_t irow = 0; irow<m_table.m_hwm; ++irow){
            if (m_table.is_cleared(irow)) continue;
            if (get_block_irow(irow) == *it_block_transfer) irows_send.push_back(irow);
        }
    }
    m_mean_work_times.assign(m_nblock, 0.0);

    /*
     * Let the rank->block and block->rank mappings reflect this reallocation
     */
    m_rank_to_blocks[irank_send].erase(it_block_transfer);
    m_rank_to_blocks[irank_recv].push_back(*it_block_transfer);
    m_block_to_rank[*it_block_transfer] = irank_recv;
    MPI_ASSERT_ALL(!m_rank_to_blocks[irank_send].empty(), "All blocks removed from rank");

    for (auto dep : m_dependents) dep->before_block_transfer(irows_send, irank_send, irank_recv);
    m_table.transfer_rows(irows_send, irank_send, irank_recv, m_recv_callbacks);
    for (auto dep : m_dependents) dep->after_block_transfer();
    MPI_ASSERT_ALL(consistent(), "block->rank map should be consistent with rank->block map");
}
