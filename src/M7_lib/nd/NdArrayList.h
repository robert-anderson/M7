//
// Created by Robert J. Anderson on 09/11/2020.
//

#ifndef M7_NDARRAYLIST_H
#define M7_NDARRAYLIST_H

#include "NdAccessor.h"

template<typename T, uint_t nind>
struct NdArrayList {
    NdFormat<nind> m_format;
    v_t<T> m_data;
    uint_t m_size;

    template<typename ...Args>
    NdArrayList(Args... shape): m_format(shape...), m_size(0){}

    void resize(uint_t size){
        m_size = size;
        m_data.resize(m_size*m_format.nelement());
    }

    void expand(uint_t nelement=1){
        resize(m_size+nelement);
    }

    NdAccessor<T, nind> operator[](const uint_t &ielement){
        ASSERT(ielement<m_size);
        return NdAccessor<T, nind>(m_data.data()+ielement*m_format.nelement(), m_format);
    }

    NdAccessor<T, nind> back(){
        return (*this)[m_size-1];
    }
};

#endif //M7_NDARRAYLIST_H
