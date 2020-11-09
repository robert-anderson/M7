//
// Created by rja on 09/11/2020.
//

#ifndef M7_NDARRAY_H
#define M7_NDARRAY_H

#include "NdAccessor.h"

template<typename T, size_t nind>
struct NdArrayBase {
    NdFormat<nind> m_format;
    std::vector<T> m_data;

    template<typename ...Args>
    NdArrayBase(Args... shape): m_format(shape...), m_data(m_format.nelement()) {}
};

template<typename T, size_t nind>
class NdArray : public NdArrayBase<T, nind>, public NdAccessor<T, nind> {
public:
    template<typename ...Args>
    NdArray(Args... shape):
            NdArrayBase<T, nind>(shape...),
            NdAccessor<T, nind>(NdArrayBase<T, nind>::m_data.data(), NdArrayBase<T, nind>::m_format) {}
};


#endif //M7_NDARRAY_H
