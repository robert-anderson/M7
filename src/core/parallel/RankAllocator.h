//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#include "src/core/field/Field.h"
#include "src/core/table/Table.h"
#include "src/core/parallel/Gatherable.h"
#include "MPIWrapper.h"
#include "Epoch.h"
#include <forward_list>
#include <algorithm>

template<typename field_t, typename hash_fn=typename field_t::hash_fn>
class RankAllocator {
    typedef typename field_t::view_t view_t;
    const size_t m_nblock;
    const size_t m_period;
    defs::inds m_block_to_rank;
    std::vector<std::forward_list<size_t>> m_rank_to_blocks;
    Gatherable<double> m_mean_wait_times;

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
    Epoch *m_vary_shift = nullptr;
    RankAllocator(const size_t& nblock, const size_t& period) :
    m_nblock(nblock), m_period(period),
    m_block_to_rank(nblock, 0ul), m_rank_to_blocks(mpi::nrank())
    {
        size_t iblock = 0ul;
        for (auto &rank:m_block_to_rank) {
            rank = iblock%mpi::nrank();
            m_rank_to_blocks[rank].push_front(iblock);
            iblock++;
        }
    }

    void update(const size_t& icycle, const double& wait_time, Table& table, field_t& field){
        if (!m_vary_shift || !*m_vary_shift) return;
        m_mean_wait_times+=wait_time;
        if (icycle%m_period) return;
        // else, this is a preparatory iteration
        auto times = m_mean_wait_times.mpi_gather();
        // this rank has done the least waiting, so it should send a block
        size_t sender = std::distance(times.begin(), std::min_element(times.begin(), times.end()));
        // this rank has done the most waiting, so it should receive a block
        size_t recver = std::distance(times.begin(), std::max_element(times.begin(), times.end()));
        if (sender==recver) return;

        if (m_rank_to_blocks[sender].empty()) mpi::stop_all("Sending rank has no blocks to send!");
        /*
         * sender will send all rows in the front block, which will become the front
         * block of the recver so all ranks must update in response to this
         */
        auto iblock_send = m_rank_to_blocks[sender].front();
        m_rank_to_blocks[sender].pop_front();
        m_rank_to_blocks[recver].push_front(iblock_send);

        // prepare vector of row indices to send
        defs::inds irows_send;
        if (mpi::i_am(sender)){
            for (size_t irow = 0; irow<table.m_hwm; ++irow){
                if (table.is_cleared(irow)) continue;
                auto key = field(irow);
                ASSERT(get_rank(key)==sender);
                if (get_block(key)==iblock_send) irows_send.push_back(irow);
            }
        }

        if (mpi::i_am(sender)) table.send_rows(irows_send, recver, m_send_callbacks);
        if (mpi::i_am(recver)) table.recv_rows(sender, m_recv_callbacks);
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
