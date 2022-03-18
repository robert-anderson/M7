//
// Created by rja on 30/03/2021.
//

#ifndef M7_NDFORMATD_H
#define M7_NDFORMATD_H

#include <utility>

#include "defs.h"
#include "NdFormat.h"

/**
 * N-dimensional formats are defined by ndim "extents"
 * these parameters can be either compile-time constant (ctc) or
 * run-time constant (rtc). The possibilities are tabulated here:
 *      ndim        extents
 *      ctc           ctc        Not implemented
 *      ctc           rtc        "nd/static"
 *      rtc           rtc        "nd/dynamic"
 */

struct NdFormatD {

    const defs::inds m_shape;
    const defs::inds m_strides;
    const std::vector<std::string> m_dim_names;
    const size_t m_nelement;
    const size_t m_nind;

    static defs::inds make_strides(const defs::inds& shape){
        const auto nind = shape.size();
        if (shape.empty()) return {};
        defs::inds strides(shape.size());
        strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            strides[nind - i] = strides[nind - i + 1] * shape[nind - i + 1];
        }
        return strides;
    }

    NdFormatD(const defs::inds& shape, std::vector<std::string> dim_names): m_shape(shape),
        m_strides(make_strides(shape)), m_dim_names(std::move(dim_names)),
        m_nelement(shape.empty() ? 1:m_strides.front()*m_shape.front()), m_nind(shape.size()){}

    NdFormatD(const defs::inds& shape): NdFormatD(shape, std::vector<std::string>(shape.size(), "")){}

    NdFormatD(size_t nind, size_t extent): NdFormatD(defs::inds(nind, extent)){}

    size_t flatten(const defs::inds& inds) const {
        size_t iflat = 0ul;
        for (size_t i=0ul; i<std::min(inds.size(), m_nind); ++i) {
            ASSERT(inds[i] < m_shape[i]);
            iflat += inds[i] * m_strides[i];
        }
        return iflat;
    }

    template<typename T>
    void decode_flat(const size_t& iflat, std::vector<T>& inds) const {
        static_assert(std::is_integral<T>::value, "index type must be integral");
        inds.clear();
        size_t remainder = iflat;
        for (size_t i=0ul; i!=m_nind; ++i){
            inds.push_back(remainder/m_strides[i]);
            remainder-=inds.back()*m_strides[i];
        }
    }


};


#endif //M7_NDFORMATD_H
