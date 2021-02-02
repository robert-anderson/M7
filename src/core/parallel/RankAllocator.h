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
    static constexpr double c_default_thresh_ratio = 0.95;
    typedef typename field_t::view_t view_t;
    Table& m_table;
    field_t& m_field;
public:
    const size_t m_nblock;
    const size_t m_period;
private:
    defs::inds m_block_to_rank;
    /*
     * push inward-transferred blocks indices to front
     * pop outward-transferred block indices to back
     */
    std::vector<std::deque<size_t>> m_rank_to_blocks;
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
        m_mean_work_times+=work_time;
        if (!icycle || icycle%m_period) return;
        // else, this is a preparatory iteration
        auto times = m_mean_work_times.mpi_gather();


        size_t irank_lazy = std::distance(times.begin(), std::min_element(times.begin(), times.end()));
        size_t irank_busy = std::distance(times.begin(), std::max_element(times.begin(), times.end()));

        // this rank has done the least work, so it should receive a block
        size_t irank_recv = irank_lazy;
        // this rank has done the most work, so it should send a block
        size_t irank_send = (irank_recv+1)%mpi::nrank();//std::distance(times.begin(), std::max_element(times.begin(), times.end()));
        // this rank has done the least work, so it should receive a block
        //size_t irank_recv = (irank_send+1)%mpi::nrank();//std::distance(times.begin(), std::min_element(times.begin(), times.end()));
        if (irank_send==irank_recv) return;
        if (times[irank_lazy] > m_threshold_ratio*times[irank_busy]) return;
        MPI_REQUIRE_ALL(!m_rank_to_blocks[irank_send].empty(),
                        "Sending rank has no blocks to send! Specified work_time value may be erroneous");
        /*
         * sender will send all rows in the front block, which will become the front
         * block of the recver so all ranks must update in response to this
         */
        auto iblock_transfer = m_rank_to_blocks[irank_send].back();
        log::info("Sending block {} from rank {} to {} on cycle {}", iblock_transfer, irank_send, irank_recv, icycle);

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

        /*
         * Let the rank->block and block->rank mappings reflect this reallocation
         */
        m_rank_to_blocks[irank_send].pop_back();
        m_rank_to_blocks[irank_recv].push_front(iblock_transfer);
        m_block_to_rank[iblock_transfer] = irank_recv;

        /*
        for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
            mpi::barrier();
            if (mpi::i_am(irank)) {
                std::cout << ">>>>> " << irank << std::endl;
                m_table.print_contents();
            }
            mpi::barrier();
        }
         */
        for (auto dep : m_dependents) dep->before_block_transfer(irows_send, irank_send, irank_recv);
        m_table.transfer_rows(irows_send, irank_send, irank_recv, m_recv_callbacks);
        for (auto dep : m_dependents) dep->after_block_transfer();
        MPI_ASSERT_ALL(consistent(), "block->rank map should be consistent with rank->block");
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
