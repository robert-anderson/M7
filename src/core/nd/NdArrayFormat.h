//
// Created by rja on 02/10/2020.
//

#ifndef M7_NDARRAYFORMAT_H
#define M7_NDARRAYFORMAT_H

#include <array>
#include <src/core/util/defs.h>

template <size_t nind>
class NdArrayFormat {
    std::array<size_t, nind> m_shape;
    std::array<size_t, nind> m_strides;

public:
    NdArrayFormat(){
        static_assert(!nind, "This ctor is only valid in the scalar case");
    }

    NdArrayFormat(size_t extent){
        m_shape.fill(extent);
        set_strides();
    }
    NdArrayFormat(std::array<size_t, nind> shape):m_shape(std::move(shape)){
        set_strides();
    }
    template<typename ...Args>
    NdArrayFormat(const size_t& first, const size_t& second, Args&&... shape){
        static_assert(2+sizeof...(Args)>=nind, "not enough arguments to specify array shape");
        static_assert(2+sizeof...(Args)<=nind, "too many arguments to specify array shape");
        set_shape(first, second, std::forward<Args>(shape)...);
        set_strides();
    }

    size_t nelement() const {
        if (!nind) return 1;
        return m_shape.front()*m_strides.front();
    }

    const size_t& extent(const size_t& i) const {
        return m_shape[i];
    }

    template<typename ...Args>
    size_t flat(Args... inds) const {
        static_assert(sizeof...(Args)==nind, "incorrect number of indices");
        return partial_offset<0>(inds...);
    }

    void decode_flat(const size_t& iflat, defs::inds& inds) const {
        ASSERT(inds.size()>=nind);
        size_t remainder = iflat;
        for (size_t i=0ul; i<nind; ++i){
            auto& ind = inds[i];
            ind = remainder/m_strides[i];
            remainder-=ind*m_strides[i];
        }
    }

private:

    void set_strides(){
        if (!nind) return;
        m_strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            m_strides[nind - i] = m_strides[nind - i + 1] * m_shape[nind - i + 1];
        }
    }

    void set_shape(){}
    template<typename ...Args>
    void set_shape(size_t first, Args... rest){
        m_shape[nind-sizeof...(rest)-1] = first;
        set_shape(rest...);
    }

    template<size_t nind_unspec>
    size_t partial_offset() const {return 0;}
    template<size_t nind_unspec, typename ...Args>
    size_t partial_offset(size_t first, Args... rest) const{
        return first*m_strides[nind-nind_unspec-1-sizeof...(Args)]+partial_offset<nind_unspec>(rest...);
    }
};


#endif //M7_NDARRAYFORMAT_H
