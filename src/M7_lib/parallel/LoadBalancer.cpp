//
// Created by rja on 17/07/22.
//

#include "LoadBalancer.h"

uint_t RankReallocator::iblock_max_le(uint_t irank_busiest, double limit) {
    uint_t iblock_best = ~0ul;
    for (auto iblock: m_rank_iblocks[irank_busiest]) {
        if ((m_block_work_figs[iblock] <= limit) && (iblock_best==~0ul ||
                                                     (m_block_work_figs[iblock] > m_block_work_figs[iblock_best]))) iblock_best = iblock;
    }
    return iblock_best;
}

bool RankReallocator::update(uint_t& iblock, uint_t& irank_dst) {
    auto it_busiest = std::max_element(m_total_work_figs.begin(), m_total_work_figs.end());
    const uint_t irank_busiest = std::distance(m_total_work_figs.begin(), it_busiest);
    DEBUG_ASSERT_GE(*it_busiest, m_perfect_work_fig,
                    "lazier than average rank should never be identified as busiest");
    auto it_laziest = std::min_element(m_total_work_figs.begin(), m_total_work_figs.end());
    irank_dst = std::distance(m_total_work_figs.begin(), it_laziest);
    DEBUG_ASSERT_LE(*it_laziest, m_perfect_work_fig,
                    "busier than average rank should never be identified as laziest");
    iblock = iblock_max_le(irank_busiest, (*it_busiest-*it_laziest)/2.0);
    if (iblock<m_block_iranks.size()) {
        DEBUG_ASSERT_EQ(m_block_iranks[iblock], irank_busiest, "chosen block does not belong to busiest rank");
        const auto diff = m_block_work_figs[iblock];
        *it_busiest -= diff;
        *it_laziest += diff;
        return true;
    }
    return false;
}

RankReallocator::RankReallocator(const uintv_t& block_iranks, const v_t<double>& block_work_figs, uint_t nrank) :
        m_block_iranks(block_iranks), m_block_work_figs(block_work_figs), m_rank_iblocks(nrank),
        m_perfect_work_fig(std::accumulate(m_block_work_figs.cbegin(), m_block_work_figs.cend(), 0.0)/nrank) {
    for (uint_t iblock=0ul; iblock<block_iranks.size(); ++iblock)
        m_rank_iblocks[block_iranks[iblock]].push_back(iblock);
    m_total_work_figs.reserve(nrank);
    for (auto& iblocks : m_rank_iblocks){
        double total_work_fig = 0.0;
        for (auto iblock: iblocks) total_work_fig+=m_block_work_figs[iblock];
        m_total_work_figs.push_back(total_work_fig);
    }
    uint_t iblock, irank_dst;
    while (update(iblock, irank_dst)) {
        m_moving_iblocks.push_back(iblock);
        m_dst_iranks.push_back(irank_dst);
    }
}
