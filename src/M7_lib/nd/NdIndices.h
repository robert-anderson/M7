//
// Created by Robert J. Anderson on 18/08/2021.
//

#ifndef M7_NDINDICES_H
#define M7_NDINDICES_H

#include "NdFormat.h"

/**
 * extension of the NdFormat which contains an enumeration of all index arrays
 * @tparam nind
 */
template<uint_t nind>
class NdIndices : NdFormat<nind> {
    using NdFormat<nind>::m_nelement;
    using NdFormat<nind>::m_strides;

    uinta_t<nind> invert(uint_t iflat) const {
        DEBUG_ASSERT_LT(iflat, m_nelement, "flat index OOB");
        if (!nind) return {};
        uinta_t<nind> inds{};
        // avoiding the usual i<nind termination condition due to unsigned type errors on some compilers
        for (uint_t i=0ul; i!=nind; ++i) {
            inds[i] = iflat / m_strides[i];
            iflat -= inds[i]*m_strides[i];
        }
        DEBUG_ASSERT_EQ(iflat, 0ul, "invert should not leave a remainder");
        return inds;
    }
    v_t<uinta_t<nind>> make_inds() const {
        v_t<uinta_t<nind>> inds;
        inds.reserve(m_nelement);
        for (uint_t iflat=0ul; iflat<m_nelement; ++iflat) inds.push_back(invert(iflat));
        return inds;
    }
public:
    const v_t<uinta_t<nind>> m_inds;

    NdIndices(): NdFormat<nind>(), m_inds(make_inds()){}

    NdIndices(const uinta_t<nind>& shape, const std::array<str_t, nind>& dim_names):
            NdFormat<nind>(shape, dim_names), m_inds(make_inds()){}

    NdIndices(const uinta_t<nind>& shape): NdFormat<nind>(shape), m_inds(make_inds()){}

    NdIndices(uint_t extent): NdFormat<nind>(extent), m_inds(make_inds()){}

    NdIndices(const NdIndices<nind>& other) : NdFormat<nind>(other), m_inds(other.m_inds){}

    uint_t size() const {
        return m_inds.size();
    }

    const uinta_t<nind>& operator[](const uint_t& i) {
        DEBUG_ASSERT_LT(i, size(), "index OOB");
        return m_inds[i];
    }
};


#endif //M7_NDINDICES_H
