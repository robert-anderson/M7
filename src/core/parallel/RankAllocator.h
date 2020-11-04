//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#if 0

#include "src/core/table/Element.h"
#include "src/core/parallel/Gatherable.h"
#include "MPIWrapper.h"
#include "Epoch.h"
#include <forward_list>
#include <algorithm>

template<typename T>
class RankAllocator {
    static_assert(std::is_base_of<Element, T>::value, "Rank allocation requires an Element-derived type");
    Epoch *m_vary_shift = nullptr;
    const size_t m_nblock;
    const size_t m_period;
    defs::inds m_block_to_rank;
    std::vector<std::forward_list<size_t>> m_rank_to_blocks;
    Gatherable<double> m_mean_wait_times;
public:
    RankAllocator(const size_t& nblock, const size_t& period, Epoch *m_vary_shift= nullptr) :
    m_vary_shift(m_vary_shift), m_nblock(nblock), m_period(period),
    m_block_to_rank(nblock, 0ul), m_rank_to_blocks(mpi::nrank())
    {
        size_t iblock = 0ul;
        for (auto &rank:m_block_to_rank) {
            rank = iblock%mpi::nrank();
            m_rank_to_blocks[rank].push_front(iblock);
            iblock++;
        }
    }

    void update(const size_t& icycle, const double& wait_time){
        if (!m_vary_shift || !*m_vary_shift) return;
        m_mean_wait_times+=wait_time;
        if (icycle%m_period) return;
        // else, this is a preparatory iteration
        auto times = m_mean_wait_times.mpi_gather();
        auto min = std::min_element(times.begin(), times.end());
        auto max = std::max_element(times.begin(), times.end());
        if (times.begin()+mpi::irank()==min){
            // this rank has done the least waiting, so it should send a block


        }
        else if (times.begin()+mpi::irank()==max){
            // this rank has done the most waiting, so it should receive a block
        }
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

    size_t get_rank(const T& key) const{
        return m_block_to_rank[key.hash()%m_nblock];
    }
};

#endif //M7_RANKALLOCATOR_H
#endif //M7_RANKALLOCATOR_H
