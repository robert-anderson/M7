//
// Created by rja on 02/10/2020.
//

#ifndef M7_NDFORMAT_H
#define M7_NDFORMAT_H

#include "src/defs.h"
#include "array"

template <size_t nind>
class NdFormat {
    std::array<size_t, nind> m_shape;
    std::array<size_t, nind> m_strides;
    size_t m_nelement = ~0ul;

    void set_nelement() {
        if (!nind) m_nelement = 1;
        else m_nelement = m_shape.front()*m_strides.front();
    }

public:
    NdFormat(){
        static_assert(!nind, "This ctor is only valid in the scalar case");
    }

    NdFormat(size_t extent){
        m_shape.fill(extent);
        set_strides();
        set_nelement();
    }
    NdFormat(std::array<size_t, nind> shape):m_shape(std::move(shape)){
        set_strides();
        set_nelement();
    }
    template<typename ...Args>
    NdFormat(const size_t& first, const size_t& second, Args&&... shape){
        static_assert(2+sizeof...(Args)>=nind, "not enough arguments to specify array shape");
        static_assert(2+sizeof...(Args)<=nind, "too many arguments to specify array shape");
        set_shape(first, second, std::forward<Args>(shape)...);
        set_strides();
        set_nelement();
    }

    NdFormat(const NdFormat<nind>& other) : NdFormat(other.shape()){}

    bool operator==(const NdFormat<nind> &other) const{
        for (size_t iind=0ul; iind<nind; ++iind){
            if (other.extent(iind)!=extent(iind)) return false;
        }
        return true;
    }

    const size_t& nelement() const {
        return m_nelement;
    }

    const size_t& extent(const size_t& i) const {
        return m_shape[i];
    }

    const size_t& stride(const size_t& i) const {
        return m_strides[i];
    }

    const std::array<size_t, nind>& shape() const {
        return m_shape;
    }

    size_t flatten(std::array<size_t, nind> inds) const {
        size_t iflat = 0ul;
        for (size_t i=0ul; i < nind; ++i) iflat+= inds[i] * m_strides[i];
        return iflat;
    }

    template<typename ...Args>
    size_t flatten(Args... inds) const {
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
public:


    /*                  nind
     *            |-------------|
     *   shape:    X  X  X  X  X
     *
     *   inds:     X  X  /  /  /
     *            |----|
     *            nind spec
     *
     */
    template<size_t nind_spec, typename ...Args>
    size_t partial_offset(size_t first, Args... rest) const{
        static_assert(1+sizeof...(rest)+nind_spec<=nind, "Indices are over-specified");
        assert(first<m_shape[nind_spec]);
        return first*m_strides[nind_spec]+partial_offset<nind_spec+1>(rest...);
    }
};


#endif //M7_NDFORMAT_H