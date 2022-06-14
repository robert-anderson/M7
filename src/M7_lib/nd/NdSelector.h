//
// Created by Robert J. Anderson on 14/02/2021.
//

#ifndef M7_NDSELECTOR_H
#define M7_NDSELECTOR_H

#include <cstring>
#include "NdFormat.h"

struct Selector {
    const size_t m_nelement_max;
    size_t m_nelement = m_nelement_max;
    size_t m_ielement = 0ul;
    Selector(size_t nelement_max): m_nelement_max(nelement_max){}
};

template<size_t nind>
struct NdSelector : Selector {
    NdFormat<nind> m_format;

    NdSelector(std::array<size_t, nind> shape):
            Selector(NdFormat<nind>(shape).nelement()), m_format(shape){}

    template<typename T, size_t nind_other>
    void copy(const T* src, const NdSelector<nind_other>& other, T* dst, size_t ndword_per_element=1){
        assert(other.m_nelement==m_nelement);
        std::memcpy(dst, src, sizeof(T)*ndword_per_element*m_nelement);
    }

    void select() {
        m_ielement = 0;
        m_nelement = m_nelement_max;
    }
    template<typename ...Args>
    void select(size_t first, Args... inds){
        static_assert(sizeof...(inds)+1 <= nind, "Invalid number of indices");
        m_ielement = m_format.template partial_offset<0ul>(first, inds...);
        m_nelement = m_format.stride(sizeof...(inds));
    }
};


#endif //M7_NDSELECTOR_H
