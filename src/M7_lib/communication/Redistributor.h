//
// Created by anderson on 27/07/2022.
//

#ifndef M7_REDISTRIBUTOR_H
#define M7_REDISTRIBUTOR_H

#include <set>
#include <numeric>
#include "M7_lib/parallel/MPIAssert.h"

/**
 * decide which block transfers to make in order to equalize load across all ranks
 */
class Redistributor {
    /**
     * map from block index to rank index
     */
    const uintv_t& m_block_iranks;
    /**
     * workload figure for each block
     */
    const v_t<double>& m_block_work_figs;
    /**
     * map from rank index to block indices
     */
    v_t<std::set<uint_t>> m_rank_iblocks;
    /**
     * total workload figure for each rank
     */
    v_t<double> m_total_work_figs;
    /**
     * total workload figure per rank with perfect load balancing
     */
    const double m_perfect_work_fig;
    /**
     * flag set to true if the rank is busier than or equal to the average load (for checks)
     */
    const v_t<bool> m_busier_or_avg;

    v_t<bool> make_busier_or_avg() const;

    /**
     * @param irank_busiest
     *  index of busiest rank (after already decided transfers)
     * @param limit
     *  maximum value of workload figure
     * @return
     *  index of block with the largest workload figure less than limit
     */
    uint_t iblock_max_lt(uint_t irank_busiest, double limit);

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

    struct Move {
        uint_t m_iblock, m_dst_irank;
        bool operator==(const Move& other) const {
            return m_iblock==other.m_iblock && m_dst_irank==other.m_dst_irank;
        }
    };
    v_t<Move> m_moves;

    Redistributor(const uintv_t& block_iranks, const v_t<double>& block_work_figs, uint_t nrank);
};


#endif //M7_REDISTRIBUTOR_H
