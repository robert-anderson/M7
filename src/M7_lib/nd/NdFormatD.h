//
// Created by Robert J. Anderson on 30/03/2021.
//

#ifndef M7_NDFORMATD_H
#define M7_NDFORMATD_H

#include <utility>

#include <M7_lib/defs.h>
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

    const uintv_t m_shape;
    const uintv_t m_strides;
    const strv_t m_dim_names;
    const uint_t m_nelement;
    const uint_t m_nind;

    static uintv_t make_strides(const uintv_t& shape){
        const auto nind = shape.size();
        if (shape.empty()) return {};
        uintv_t strides(shape.size());
        strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            strides[nind - i] = strides[nind - i + 1] * shape[nind - i + 1];
        }
        return strides;
    }

    NdFormatD(const uintv_t& shape, strv_t dim_names): m_shape(shape),
                                                                               m_strides(make_strides(shape)), m_dim_names(std::move(dim_names)),
                                                                               m_nelement(shape.empty() ? 1ul: m_strides.front()*m_shape.front()), m_nind(shape.size()){}

    NdFormatD(const uintv_t& shape): NdFormatD(shape, strv_t(shape.size(), "")){}

    NdFormatD(uint_t nind, uint_t extent): NdFormatD(uintv_t(nind, extent)){}

    uint_t flatten(const uintv_t& inds) const {
        uint_t iflat = 0ul;
        for (uint_t i=0ul; i<std::min(inds.size(), m_nind); ++i) {
            ASSERT(inds[i] < m_shape[i]);
            iflat += inds[i] * m_strides[i];
        }
        return iflat;
    }

    template<typename T>
    void decode_flat(const uint_t& iflat, std::vector<T>& inds) const {
        static_assert(std::is_integral<T>::value, "index type must be integral");
        inds.clear();
        uint_t remainder = iflat;
        for (uint_t i = 0ul; i != m_nind; ++i) {
            inds.push_back(remainder / m_strides[i]);
            remainder -= inds.back() * m_strides[i];
        }
    }
};


struct NdEnumerationD : NdFormatD {
private:
    const std::vector<uintv_t> m_inds;

    static std::vector<uintv_t> make_inds(const NdFormatD& format) {
        using namespace basic_foreach::rtnd;
        std::vector<uintv_t> out(format.m_nelement);
        uint_t i=0ul;
        auto fn = [&out, &i](const uintv_t& inds){
            out[i] = inds;
            ++i;
        };
        Unrestricted(format.m_shape).loop(fn);
        DEBUG_ASSERT_EQ(i, format.m_nelement, "not all index arrays generated");
        return out;
    }

public:
    NdEnumerationD(const NdFormatD& format): NdFormatD(format), m_inds(make_inds(format)){}

    const uintv_t& operator[](uint_t i) const {
        DEBUG_ASSERT_LT(i, m_inds.size(), "index OOB");
        return m_inds[i];
    }

};


#endif //M7_NDFORMATD_H
