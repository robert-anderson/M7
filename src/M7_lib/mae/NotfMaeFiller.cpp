//
// Created by anderson on 09/08/2023.
//

#include "NotfMaeFiller.h"

v_t<uintv_t> NotfMaeFiller::make_occ_bitsets(uint_t displ, uint_t count) const {
    v_t<uintv_t> bitsets;
    const auto nspinorb = m_hist.m_row.m_mbf.m_basis.m_nspinorb;
    for (uint_t ispinorb=0ul; ispinorb < nspinorb; ++ispinorb) {
        uintv_t inds;
        auto row = m_hist.m_row;
        for (row.restart(displ); row.in_range(displ + count); ++row) {
            if (row.m_mbf.get(ispinorb)) inds.push_back(row.index());
        }
        bitsets.emplace_back(bitset_isect::make_bitset(inds));
    }
    return bitsets;
}

v_t<uintv_t> NotfMaeFiller::make_occ_bitsets() const {
    return make_occ_bitsets(0, m_hist.nrow_in_use());
}
