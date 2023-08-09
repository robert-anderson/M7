//
// Created by anderson on 09/08/2023.
//

#include "NotfMaeFiller.h"

// constructor 1?
v_t<uintv_t> NotfMaeFiller::make_occ_bitsets(uint_t displ, uint_t count) const {
    // initialised to empty
    v_t<uintv_t> bitsets;  // {} ?
    const auto nspinorb = m_hist.m_row.m_mbf.m_basis.m_nspinorb;  // get a spinorb
    for (uint_t ispinorb=0ul; ispinorb < nspinorb; ++ispinorb) {  // ul = unsigned long
        uintv_t inds;
        auto row = m_hist.m_row;
        // get chunk of histdets on rank
        for (row.restart(displ); row.in_range(displ + count); ++row) {
            // if spinorb occ append to set
            if (row.m_mbf.get(ispinorb)) inds.push_back(row.index());
        }
        // append vector to vector of vectors
        bitsets.emplace_back(bitset_isect::make_bitset(inds));
    }
    return bitsets;
}

// constructor 2?
v_t<uintv_t> NotfMaeFiller::make_occ_bitsets() const {
    return make_occ_bitsets(0, m_hist.nrow_in_use());
}
