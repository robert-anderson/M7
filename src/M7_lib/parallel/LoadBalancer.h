//
// Created by rja on 17/07/22.
//

#ifndef M7_LOADBALANCER_H
#define M7_LOADBALANCER_H

#include <numeric>
#include <algorithm>
#include "MPIAssert.h"

/**
 * decide which block transfers to make in order to equalize load across all ranks
 */
class RankReallocator {
    /**
     * map from block index to rank index
     */
    const uintv_t& m_block_iranks;
    /**
     * workload figure for each block
     */
    const v_t<double>& m_block_work_figs;
    /**
     * vectors of block indices for each rank
     */
    v_t<uintv_t> m_rank_iblocks;
    /**
     * total workload figure for each rank
     */
    v_t<double> m_total_work_figs;
    /**
     * total workload figure per rank with perfect load balancing
     */
    const double m_perfect_work_fig;

    /**
     * @param irank_busiest
     *  index of busiest rank (after already decided transfers)
     * @param limit
     *  maximum value of workload figure
     * @return
     *  index of block with the largest workload figure less than limit
     */
    uint_t iblock_max_le(uint_t irank_busiest, double limit);

    /**
     * do one iteration of the reallocation algorithm
     * @param iblock
     *  index of the block to transfer
     * @param irank_dst
     *  index of MPI rank to which block iblock is destined
     * @return
     *  false if no further transfers could be found
     */
    bool update(uint_t& iblock, uint_t& irank_dst);
public:

    /**
     * block indices to move
     */
    uintv_t m_moving_iblocks;
    /**
     * indices of lazy ranks for which each moving block is destined
     */
    uintv_t m_dst_iranks;

    RankReallocator(const uintv_t& block_iranks, const v_t<double>& block_work_figs, uint_t nrank);
};


struct LoadBalancerBase {
    /**
     * workload figures for each block
     */
    v_t<double> m_block_work_figs;

    uintv_t m_block_to_rank;

    void accumulate_work_fig(uint_t iblock, double cost){
        m_block_work_figs[iblock]+=cost;
    }

};


struct LoadBalancer {

};


#endif //M7_LOADBALANCER_H
