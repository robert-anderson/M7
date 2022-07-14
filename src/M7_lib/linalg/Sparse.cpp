//
// Created by Robert J. Anderson on 14/06/2020.
//

#include "Sparse.h"
#include "M7_lib/util/String.h"

uintv_t sparse::fixed::Base::make_counts(const sparse::dynamic::Base &src) {
    uintv_t counts;
    counts.reserve(src.nrow());
    for (uint_t irow=0ul; irow<src.nrow(); ++irow) counts.push_back(src.nentry(irow));
    return counts;
}

uintv_t sparse::fixed::Base::make_displs(const uintv_t &counts) {
    uintv_t displs;
    /*
     * include the last one so that the cend iterator for the last entry is directly accessible
     */
    displs.reserve(counts.size()+1);
    displs.push_back(0ul);
    for (auto it = counts.cbegin(); it!=counts.cend(); ++it) displs.push_back(displs.back()+*it);
    return displs;
}

sparse::fixed::Base::Base(const uintv_t &counts) :
    m_nrow(counts.size()), m_max_nentry(counts.empty() ? 0ul : *std::max(counts.cbegin(), counts.cend())),
    m_displs(make_displs(counts)), m_nentry(m_displs.back()) {}

sparse::fixed::Base::Base(const sparse::dynamic::Base &src) : Base(make_counts(src)){
    DEBUG_ASSERT_EQ(m_nentry, src.nentry(),
                    "number of entries calculated from offsets should match that counted by the dynamic::Network instance");
    DEBUG_ASSERT_EQ(m_max_nentry, src.max_col_ind()+1,
                    "max number of entries in row calculated from offsets should match dynamic::Network::m_icol_max+1");
}