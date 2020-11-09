//
// Created by rja on 09/11/2020.
//

#ifndef M7_NDARRAYLIST_H
#define M7_NDARRAYLIST_H

#include "NdAccessor.h"

template<typename T, size_t nind>
struct NdArrayList {
    NdFormat<nind> m_format;
    std::vector<T> m_data;
    size_t m_size;

    template<typename ...Args>
    NdArrayList(Args... shape): m_format(shape...), m_size(0){}

    void resize(size_t size){
        m_size = size;
        m_data.resize(m_size*m_format.nelement());
    }

    NdAccessor<T, nind> operator[](const size_t &ielement){
        ASSERT(ielement<m_size);
        return NdAccessor<T, nind>(m_data.data()+ielement*m_format.nelement(), m_format);
    }
};

#endif //M7_NDARRAYLIST_H
