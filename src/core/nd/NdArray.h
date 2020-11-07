//
// Created by rja on 12/10/2020.
//

#ifndef M7_NDARRAY_H
#define M7_NDARRAY_H

#include "NdFormat.h"

template <typename T, size_t nind>
struct NdArray {
    NdFormat<nind> m_format;
    std::vector<T> m_data;

    template<typename ...Args>
    NdArray(Args... shape): m_format(shape...), m_data(m_format.nelement()){}

    template<typename ...Args>
    T& operator()(Args... inds){
        return m_data[m_format.flatten(inds...)];
    }

    template<typename ...Args>
    const T& operator()(Args... inds) const {
        return m_data[m_format.flatten(inds...)];
    }

    size_t nelement() const {
        return m_format.nelement();
    }

    void zero() {
        m_data.assign(m_data.size(), T{});
    }
};


#endif //M7_NDARRAY_H
