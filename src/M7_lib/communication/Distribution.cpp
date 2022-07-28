//
// Created by anderson on 27/07/2022.
//

#include "Distribution.h"

Distribution::Distribution(size_t nblock) {
    REQUIRE_GE(nblock, mpi::nrank(), "number of blocks may not be less than the number of ranks");
    m_block_iranks.reserve(nblock);
    /*
     * initialize distribution evenly
     */
    for (uint_t irank=0ul; irank < mpi::nrank(); ++irank){
        m_nblock_local = mpi::evenly_shared_count(nblock);
        for (uint_t iblock = 0ul; iblock < m_nblock_local; ++iblock) {
            m_block_iranks.push_back(irank);
        }
    }
    DEBUG_ASSERT_EQ(m_block_iranks.size(), nblock, "error in initial block allocation");
}

void Distribution::update(const Redistributor &redist) {
    for (auto move: redist.m_moves) {
        auto& irank = m_block_iranks[move.m_iblock];
        if (mpi::i_am(irank)) --m_nblock_local;
        irank = move.m_dst_irank;
        if (mpi::i_am(irank)) ++m_nblock_local;
    }
}
