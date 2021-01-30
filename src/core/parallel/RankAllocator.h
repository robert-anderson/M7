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

template<typename field_t, typename hash_fn=typename field_t::hash_fn>
class RankAllocator {
    static constexpr double c_default_thresh_ratio = 0.9;
    typedef typename field_t::view_t view_t;
    Table& m_table;
    field_t& m_field;
public:
    const size_t m_nblock;
    const size_t m_period;
private:
    defs::inds m_block_to_rank;
    std::vector<std::forward_list<size_t>> m_rank_to_blocks;
    Gatherable<double> m_mean_work_times;
    /*
     * if the least productive rank is working for this proportion of the
     * time worked by the most productive rank, then move a block
     */
    double m_threshold_ratio = c_default_thresh_ratio;

public:
    struct Dynamic {
        typedef RankAllocator<field_t, hash_fn> ra_t;
        ra_t& m_ra;
        typename std::list<Dynamic*>::iterator m_it;
        Dynamic(ra_t& ra): m_ra(ra), m_it(m_ra.add_dependent(this)) {}
        ~Dynamic() {
            m_ra.m_dependents.erase(m_it);
        }

        // called before send
        virtual void on_row_send_(size_t irow) = 0;
        // called after recv
        virtual void on_row_recv_(size_t irow) = 0;
    };

private:
    std::list<Dynamic*> m_dependents;

    /*
     * for each dependent, append a lambda to a list which is then passed to the Table's
     * send and recv methods
     */
    Table::cb_list_t m_send_callbacks;
    Table::cb_list_t m_recv_callbacks;

    typename std::list<Dynamic*>::iterator add_dependent(Dynamic* dependent){
        m_dependents.push_back(dependent);
        auto it = m_dependents.end();

        // refresh callback lists
        m_send_callbacks.clear();
        m_recv_callbacks.clear();
        for (auto ptr : m_dependents) {
            m_send_callbacks.push_back([ptr](const size_t &irow) { ptr->on_row_send_(irow); });
            m_recv_callbacks.push_back([ptr](const size_t &irow) { ptr->on_row_recv_(irow); });
        }
        return --it;
    }

public:
    RankAllocator(Table& table, field_t& field, const size_t& nblock, const size_t& period) :
    m_table(table), m_field(field), m_nblock(nblock), m_period(period),
    m_block_to_rank(nblock, 0ul), m_rank_to_blocks(mpi::nrank())
    {
        size_t iblock = 0ul;
        for (auto &rank:m_block_to_rank) {
            rank = iblock%mpi::nrank();
            m_rank_to_blocks[rank].push_front(iblock);
            iblock++;
        }
    }

    void update(const size_t& icycle, const double& work_time){
        log::info_("work time: {}", work_time);
        m_mean_work_times+=work_time;
        if (!icycle || icycle%m_period) return;
        // else, this is a preparatory iteration
        auto times = m_mean_work_times.mpi_gather();
        // this rank has done the most work, so it should send a block
        size_t irank_send = std::distance(times.begin(), std::max_element(times.begin(), times.end()));
        // this rank has done the least work, so it should receive a block
        size_t irank_recv = std::distance(times.begin(), std::min_element(times.begin(), times.end()));
        if (irank_send==irank_recv) return;
        if (times[irank_recv] > m_threshold_ratio*times[irank_send]) return;
        MPI_REQUIRE_ALL(!m_rank_to_blocks[irank_send].empty(), "Sending rank has no blocks to send!");
        /*
         * sender will send all rows in the front block, which will become the front
         * block of the recver so all ranks must update in response to this
         */
        auto iblock_transfer = m_rank_to_blocks[irank_send].front();
        log::info("Sending block {} from rank {} to rank {}", iblock_transfer, irank_send, irank_recv);
        m_rank_to_blocks[irank_send].pop_front();
        m_rank_to_blocks[irank_recv].push_front(iblock_transfer);

        // prepare vector of row indices to send
        defs::inds irows_send;
        if (mpi::i_am(irank_send)){
            for (size_t irow = 0; irow<m_table.m_hwm; ++irow){
                if (m_table.is_cleared(irow)) continue;
                auto key = m_field(irow);
                ASSERT(get_rank(key)==irank_send);
                if (get_block(key) == iblock_transfer) irows_send.push_back(irow);
            }
        }

        if (mpi::i_am(irank_send)) m_table.send_rows(irows_send, irank_recv, m_send_callbacks);
        if (mpi::i_am(irank_recv)) m_table.recv_rows(irank_send, m_recv_callbacks);
        m_mean_work_times = 0.0;
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
