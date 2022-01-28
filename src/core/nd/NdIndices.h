//
// Created by rja on 18/08/2021.
//

#ifndef M7_NDINDICES_H
#define M7_NDINDICES_H

#include "NdFormat.h"

/**
 * extension of the NdFormat which contains an enumeration of all index arrays
 * @tparam nind
 */
template<size_t nind>
class NdIndices : NdFormat<nind> {
    using NdFormat<nind>::m_nelement;
    using NdFormat<nind>::m_strides;

    std::array<size_t, nind> invert(size_t iflat) const {
        DEBUG_ASSERT_LT(iflat, m_nelement, "flat index OOB");
        if (!nind) return {};
        std::array<size_t, nind> inds{};
        // avoiding the usual i<nind termination condition due to unsigned type errors on some compilers
        for (size_t i=0ul; i!=nind; ++i) {
            inds[i] = iflat / m_strides[i];
            iflat -= inds[i]*m_strides[i];
        }
        DEBUG_ASSERT_EQ(iflat, 0ul, "invert should not leave a remainder");
        return inds;
    }
    std::vector<std::array<size_t, nind>> make_inds() const {
        std::vector<std::array<size_t, nind>> inds;
        inds.reserve(m_nelement);
        for (size_t iflat=0ul; iflat<m_nelement; ++iflat) inds.push_back(invert(iflat));
        return inds;
    }
public:
    const std::vector<std::array<size_t, nind>> m_inds;

    NdIndices(): NdFormat<nind>(), m_inds(make_inds()){}

    NdIndices(const std::array<size_t, nind>& shape, const std::array<std::string, nind>& dim_names):
            NdFormat<nind>(shape, dim_names), m_inds(make_inds()){}

    NdIndices(const std::array<size_t, nind>& shape): NdFormat<nind>(shape), m_inds(make_inds()){}

    NdIndices(size_t extent): NdFormat<nind>(extent), m_inds(make_inds()){}

    NdIndices(const NdIndices<nind>& other) : NdFormat<nind>(other), m_inds(other.m_inds){}

    size_t size() const {
        return m_inds.size();
    }

    const std::array<size_t, nind>& operator[](const size_t& i) {
        DEBUG_ASSERT_LT(i, size(), "index OOB");
        return m_inds[i];
    }
};


#endif //M7_NDINDICES_H
