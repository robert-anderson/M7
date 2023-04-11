//
// Created by anderson on 27/07/2022.
//

#include "Redistributor.h"

v_t<bool> Redistributor::make_busier_or_avg() const {
    v_t<bool> flags;
    flags.reserve(m_total_work_figs.size());
    for (auto v: m_total_work_figs) flags.push_back(v>=m_perfect_work_fig);
    return flags;
}

uint_t Redistributor::iblock_max_lt(uint_t irank_busiest, double limit) {
    uint_t iblock_best = ~0ul;
    for (auto iblock: m_rank_iblocks[irank_busiest]) {
        if ((m_block_work_figs[iblock] < limit) && (iblock_best==~0ul ||
                                                    (m_block_work_figs[iblock] > m_block_work_figs[iblock_best]))) iblock_best = iblock;
    }
    return iblock_best;
}

bool Redistributor::update(uint_t &iblock, uint_t &irank_dst) {
    auto it_src = std::max_element(m_total_work_figs.begin(), m_total_work_figs.end());
    const uint_t irank_src = std::distance(m_total_work_figs.begin(), it_src);
    auto it_dst = std::min_element(m_total_work_figs.begin(), m_total_work_figs.end());
    irank_dst = std::distance(m_total_work_figs.begin(), it_dst);
    DEBUG_ASSERT_NE(irank_src, irank_dst, "rank should never transfer to itself");

    // edge case where one of the extreme ranks has exactly the perfect load
    if (*it_src == m_perfect_work_fig || *it_dst == m_perfect_work_fig) {
        DEBUG_ASSERT_TRUE(*it_src == m_perfect_work_fig || *it_dst == m_perfect_work_fig,
                          "both should be perfectly distributed");
        return false;
    }

    DEBUG_ASSERT_GT(*it_src, m_perfect_work_fig,
                    "lazier than average rank should never be identified as busiest");
    DEBUG_ASSERT_LT(*it_dst, m_perfect_work_fig,
                    "busier than average rank should never be identified as laziest");

    /*
     * find the index of the block with the largest weight less than half the difference between the work figures of
     * the current busiest and laziest ranks. this will be the most efficient move. if the work figure is greater
     * than this limit, the lazy rank will become busier than average
     */
    iblock = iblock_max_lt(irank_src, (*it_src - *it_dst) / 2.0);
    if (iblock<m_block_iranks.size()) {
        DEBUG_ASSERT_FALSE(m_rank_iblocks[irank_src].empty(), "source rank should have at least one block");
        const auto diff = m_block_work_figs[iblock];
        *it_src -= diff;
        *it_dst += diff;
        m_rank_iblocks[irank_src].erase(iblock);
        m_rank_iblocks[irank_dst].insert(iblock);
        return true;
    }
    return false;
}

Redistributor::Redistributor(const uintv_t &block_iranks, const v_t<double> &block_work_figs, uint_t nrank) :
        m_block_iranks(block_iranks), m_block_work_figs(block_work_figs), m_rank_iblocks(nrank),
        m_perfect_work_fig(std::accumulate(m_block_work_figs.cbegin(), m_block_work_figs.cend(), 0.0)/nrank),
        m_busier_or_avg(make_busier_or_avg()){
    {
        /*
         * construct the inverse one-to-many map from ranks to blocks
         */
        uint_t iblock = 0ul;
        for (auto& irank: block_iranks) m_rank_iblocks[irank].insert(iblock++);
    }

    m_total_work_figs.resize(nrank);
    for (uint_t iblock=0ul; iblock<m_block_iranks.size(); ++iblock) {
        m_total_work_figs[m_block_iranks[iblock]] += m_block_work_figs[iblock];
        DEBUG_ASSERT_GE(m_block_work_figs[iblock], 0.0, "block work figure should be non-negative");
    }
    if (nrank == 1ul) return;
    uint_t iblock, irank_dst;
    while (update(iblock, irank_dst)) m_moves.push_back({iblock, irank_dst});
}
