//
// Created by rja on 02/10/2020.
//

#ifndef M7_NDFORMAT_H
#define M7_NDFORMAT_H

#include <src/core/enumerator/Enumerator.h>
#include "src/defs.h"
#include "array"
#include "algorithm"

struct NdFormatBase{
    virtual std::string to_string() const = 0;
};

template <size_t nind>
class NdFormat : NdFormatBase{
    std::array<size_t, nind> m_shape;
    std::array<size_t, nind> m_strides;
    std::array<std::string, nind> m_dim_names;
    size_t m_nelement = ~0ul;

    void set_nelement() {
        if (!nind) m_nelement = 1;
        else m_nelement = m_shape.front()*m_strides.front();
    }

public:
    NdFormat(){
        static_assert(!nind, "This ctor is only valid in the scalar case");
        m_nelement = 1ul;
    }

    NdFormat(size_t extent){
        m_shape.fill(extent);
        set_strides();
        set_nelement();
    }

    NdFormat(const std::array<size_t, nind>& shape, const std::array<std::string, nind>& dim_names={}):
        m_shape(shape), m_dim_names(dim_names){
        set_strides();
        set_nelement();
        ASSERT(m_nelement!=~0ul);
    }

    NdFormat(const NdFormat<nind>& other) : NdFormat(other.m_shape, other.m_dim_names){}

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

    defs::inds shape_vector() const {
        defs::inds tmp;
        tmp.assign(m_shape.cbegin(), m_shape.cend());
        return tmp;
    }

    const std::array<std::string, nind>& dim_names() const {
        return m_dim_names;
    }
    std::vector<std::string> dim_names_vector() const {
        std::vector<std::string> tmp;
        tmp.assign(m_dim_names.cbegin(), m_dim_names.cend());
        return tmp;
    }

    size_t flatten(std::array<size_t, nind> inds) const {
        size_t iflat = 0ul;
        for (size_t i=0ul; i < nind; ++i) {
            ASSERT(inds[i]<m_shape[i]);
            iflat+= inds[i] * m_strides[i];
        }
        return iflat;
    }

    template<typename ...Args>
    size_t flatten(Args... inds) const {
        static_assert(sizeof...(Args)==nind, "incorrect number of indices");
        return partial_offset<0>(inds...);
    }

    void decode_flat(const size_t& iflat, std::array<size_t, nind>& inds) const {
        size_t remainder = iflat;
        for (size_t i=0ul; i<nind; ++i){
            auto& ind = inds[i];
            ind = remainder/m_strides[i];
            remainder-=ind*m_strides[i];
        }
    }

    std::string to_string() const override {
        if (!nind) return "scalar";
        std::string tmp;
        auto size = nind;
        if (m_dim_names==std::array<std::string, nind>{}){
            for (size_t i=0ul; i<size; ++i) {
                tmp+=" "+std::to_string(m_shape[i]);
            }
        }
        else {
            for (size_t i=0ul; i<size; ++i) {
                tmp+=m_dim_names[i]+" ("+std::to_string(m_shape[i])+") ";
            }
        }
        return tmp;
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
        ASSERT(first<m_shape[nind_spec]);
        return first*m_strides[nind_spec]+partial_offset<nind_spec+1>(rest...);
    }

    template<size_t nind_spec>
    size_t partial_offset(std::array<size_t, nind_spec> inds) const {
        static_assert(nind_spec<=nind, "Too many indices specified");
        size_t iflat = 0ul;
        for (size_t i=0ul; i < nind_spec; ++i) {
            ASSERT(inds[i]<m_shape[i]);
            iflat+= inds[i] * m_strides[i];
        }
        return iflat;
    }
};

template <size_t nind>
struct NdEnumeration {
    const NdFormat<nind> m_format;
    const std::vector<std::array<size_t, nind>> m_inds;

    static std::vector<std::array<size_t, nind>> make_inds(const NdFormat<nind>& format) {
        std::vector<std::array<size_t, nind>> out(format.nelement());
        enums::PermutationsWithRepetition e(format.shape_vector());
        size_t iinds = 0ul;
        while (e.next()){
            for (size_t i=0; i<nind; ++i) out[iinds][i] = e[i];
            ++iinds;
        }
        ASSERT(iinds==format.nelement());
        return out;
    }

    NdEnumeration(const NdFormat<nind>& format): m_format(format), m_inds(make_inds(format)){}

    const size_t& nelement() const {
        return m_format.nelement();
    }

    const std::array<size_t, nind>& operator[](const size_t& i){
        return m_inds[i];
    }

    operator const NdFormat<nind>&() const {
        return m_format;
    }

};


#endif //M7_NDFORMAT_H