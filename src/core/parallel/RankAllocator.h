//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#include "src/core/field/Field.h"
#include "src/core/table/Table.h"
#include "src/core/parallel/Gatherable.h"
#include "src/core/io/Logging.h"
#include "MPIWrapper.h"
#include "Epoch.h"
#include <forward_list>
#include <algorithm>
#include <src/core/util/Timer.h>

template<typename field_t, typename hash_fn=typename field_t::hash_fn>
class RankAllocator {
    static constexpr double c_default_acceptable_imbalance = 0.01;
    static constexpr size_t c_nnull_updates_deactivate = 20;
    typedef typename field_t::view_t view_t;
    Table& m_table;
    field_t& m_field;
public:
    /*
     * we should be able to deactivate the load balancing once certain externally
     * specified conditions are met
     */
    bool m_active = true;
    /*
     * total number of blocks across all ranks
     */
    const size_t m_nblock;
    /*
     * number of cycles between reallocation attempts
     */
    const size_t m_period;
private:
    defs::inds m_block_to_rank;
    /*
     * push inward-transferred blocks indices to front
     * pop outward-transferred block indices to back
     */
    std::vector<std::list<size_t>> m_rank_to_blocks;
    std::vector<double> m_mean_work_times;
    std::vector<double> m_gathered_times;
    /*
     * if the least productive rank is working for this proportion of the
     * time worked by the most productive rank, then move a block
     */
    double m_acceptable_imbalance = c_default_acceptable_imbalance;
    /*
     * If the load balancing algorithm is stable, a situation should be reached in which
     * the update method finds no reason to reallocate blocks. If this happens for a
     * certain m_nnull_updates_deactivate periods, m_active is automatically set to false.
     */
    size_t m_nnull_updates_deactivate = c_nnull_updates_deactivate;
    size_t m_nnull_updates = 0ul;
    /*
     * moving blocks with a large overhead is often counterproductive to achieving good
     * load balance, so if a block has a mean_work_time above the
     * m_freeze_ratio * local_work_time, it is not to be moved.
     *
     * The default value is twice the expected time if all blocks were equally expensive
     */
    double m_freeze_ratio = 2.0/m_nblock;

public:
    struct Dynamic {
        typedef RankAllocator<field_t, hash_fn> ra_t;
        ra_t& m_ra;
        typename std::list<Dynamic*>::iterator m_it;
        Dynamic(ra_t& ra): m_ra(ra), m_it(m_ra.add_dependent(this)) {}
        ~Dynamic() {
            m_ra.erase_dependent(this);
        }

        // called before block transfer is performed
        virtual void before_block_transfer(const defs::inds& irows_send, size_t irank_send, size_t irank_recv) = 0;
        // called after row is inserted into the MappedTable on the recving rank
        virtual void on_row_recv_(size_t irow) = 0;
        // called after block transfer is performed
        virtual void after_block_transfer() = 0;
    };

private:
    std::list<Dynamic*> m_dependents;
    /*
     * for each dependent, append a lambda to a list which is then called in the update method
     */
    std::list<Table::recv_cb_t> m_recv_callbacks;

    void refresh_callback_list(){
        m_recv_callbacks.clear();
        for (auto ptr : m_dependents) {
            m_recv_callbacks.push_back([ptr](size_t irow) { ptr->on_row_recv_(irow); });
        }
    }

    typename std::list<Dynamic*>::iterator add_dependent(Dynamic* dependent){
        m_dependents.push_back(dependent);
        auto it = m_dependents.end();
        refresh_callback_list();
        return --it;
    }

    void erase_dependent(Dynamic* dependent){
        m_dependents.erase(dependent->m_it);
        refresh_callback_list();
    }

public:
    RankAllocator(Table& table, field_t& field, const size_t& nblock, const size_t& period) :
    m_table(table), m_field(field), m_nblock(nblock), m_period(period),
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

    void record_work_time(const size_t& irow, const Timer& work_time){
        m_mean_work_times[get_block(m_field(irow))]+=work_time;
    }

    size_t get_nskip_(const size_t total_local_time) const {
        /*
         * When selecting a block to transfer, we need to skip over all frozen blocks.
         *
         * only called on the sending rank since the block-wise work times are not
         * shared with all ranks.
         *
         * The selected block is the first unfrozen block from the front of the list
         */
        const auto& list = m_rank_to_blocks[mpi::irank()];
        size_t nskip=0ul;
        for (auto it=list.cbegin(); it!=list.cend(); ++it){
            if (m_mean_work_times[*it] < m_freeze_ratio*total_local_time) return nskip;
            ++nskip;
        }
        log::warn_("All blocks on sending rank are frozen! Sending expensive block, may upset load balance");
        /*
         * If we reach here, all blocks are frozen, and so we have no choice but to send one of them,
         * choose the least expensive one.
         */
        double least = std::numeric_limits<double>::max();
        nskip = 0ul;
        for (auto it=list.cbegin(); it!=list.cend(); ++it){
            if (m_mean_work_times[*it] < least) nskip = std::distance(list.cbegin(), it);
        }
        return nskip;
    }

    void update(const size_t& icycle){
        if (!m_active) return;
        if (!icycle || icycle%m_period) return;
        // else, this is a preparatory iteration
        auto local_time = std::accumulate(m_mean_work_times.begin(), m_mean_work_times.end(), 0.0);
        mpi::all_gather(local_time, m_gathered_times);

        auto it_lazy = std::min_element(m_gathered_times.begin(), m_gathered_times.end());
        auto it_busy = std::max_element(m_gathered_times.begin(), m_gathered_times.end());
        MPI_REQUIRE_ALL(*it_lazy>0.0 && *it_busy>0.0,
                        "One or more ranks appear to have done no work. record_work_time must be called by application");

        // this rank has done the least work, so it should receive a block
        size_t irank_recv = std::distance(m_gathered_times.begin(), it_lazy);
        // this rank has done the most work, so it should send a block
        size_t irank_send = std::distance(m_gathered_times.begin(), it_busy);

        if (irank_send==irank_recv) return;
        auto realloc_ratio = 1-m_acceptable_imbalance;
        if (m_gathered_times[irank_recv] > realloc_ratio * m_gathered_times[irank_send]) {
            ++m_nnull_updates;
            if(m_nnull_updates>=m_nnull_updates_deactivate) {
                log::info(
                        "Load imbalance has been below {}% for over {} periods of {} cycles. "
                        "Deactivating dynamic rank allocation.",
                        m_acceptable_imbalance * 100, m_nnull_updates, m_period);
                m_active = false;
                // rezero the counter in case of later reactivation
                m_nnull_updates = 0ul;
            }
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
        if (mpi::i_am(irank_send)) nskip = get_nskip_(local_time);
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
                auto key = m_field(irow);
                ASSERT(get_rank(key)==irank_send);
                if (get_block(key) == *it_block_transfer) irows_send.push_back(irow);
            }
        }

        /*
         * Let the rank->block and block->rank mappings reflect this reallocation
         */
        m_rank_to_blocks[irank_send].erase(it_block_transfer);
        m_rank_to_blocks[irank_recv].push_front(*it_block_transfer);
        m_block_to_rank[*it_block_transfer] = irank_recv;

        for (auto dep : m_dependents) dep->before_block_transfer(irows_send, irank_send, irank_recv);
        m_table.transfer_rows(irows_send, irank_send, irank_recv, m_recv_callbacks);
        for (auto dep : m_dependents) dep->after_block_transfer();
        MPI_ASSERT_ALL(consistent(), "block->rank map should be consistent with rank->block map");
        m_mean_work_times.assign(m_nblock, 0.0);
    }

    /**
     * checks that the m_block_to_rank map is consistent with its inverse: m_rank_to_blocks
     * @return true if maps are consistent
     */
    bool consistent(){
        for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
            for (auto& iblock: m_rank_to_blocks[irank]){
                if (m_block_to_rank[iblock]!=irank) return false;
            }
        }
        return true;
    }

    inline size_t get_block(const view_t& key) const{
        return hash_fn()(key)%m_nblock;
    }

    inline size_t get_rank(const view_t& key) const{
        return m_block_to_rank[get_block(key)];
    }
};

namespace ra {
    using Onv = RankAllocator<fields::Onv<>>;
}


#endif //M7_RANKALLOCATOR_H
