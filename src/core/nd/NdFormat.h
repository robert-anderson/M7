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

/**
 * Represents a multidimensional indexing scheme where the number of dimensions is known at compile time
 * @tparam nind
 *  number of dimensions.
 */
template <size_t nind>
class NdFormat : NdFormatBase{
    /**
     * the extents of each dimension
     */
    std::array<size_t, nind> m_shape;
    /**
     * elements of the stride are computed as products of elements of the shape. this allows for the conversion of an
     * array containing multidimensional indices to a single integer ("flat" index) via a scalar product
     */
    std::array<size_t, nind> m_strides;
    /**
     * label each dimension
     */
    std::array<std::string, nind> m_dim_names;
    /**
     * the number of elements in the indexing scheme. the largest value a flat index can take in this scheme is one less
     * than m_nelement
     */
    size_t m_nelement = ~0ul;

    /**
     * helper function to store the number of elements
     */
    void set_nelement() {
        if (!nind) m_nelement = 1;
        else m_nelement = m_shape.front()*m_strides.front();
    }

public:
    NdFormat(){
        static_assert(!nind, "This ctor is only valid in the scalar case");
        m_nelement = 1ul;
    }

    /**
     * construct a scheme where all extents are the same
     * @param extent
     *  value to copy to all elements of the shape
     */
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

    /**
     * get another NdFormat where the nind-nind_remain most minor indices have been discarded. Individual elements of
     * the result index subarrays of this object
     * @tparam nind_remain
     *  number of major (right-most) indices to keep
     * @return
     *  new NdFormat object where the major part of the shape is retained
     */
    template <size_t nind_remain = nind-1>
    typename std::enable_if<(nind_remain<nind), NdFormat<nind_remain>>::type
    major_dims() const {
        constexpr auto nind_diff = nind-nind_remain;
        std::array<size_t, nind_remain> shape;
        std::copy(m_shape.cbegin(), m_shape.cend()-nind_diff, shape.begin());
        std::array<std::string, nind_remain> dim_names;
        std::copy(m_dim_names.cbegin(), m_dim_names.cend()-nind_diff, dim_names.begin());
        return {shape, dim_names};
    }
    /**
     * same as major_dims except the discarded indices are major
     * @tparam nind_remain
     *  number of minor (left-most) indices to keep
     * @return
     *  new NdFormat object where the minor part of the shape is retained
     */
    template <size_t nind_remain = nind-1>
    typename std::enable_if<(nind_remain<nind), NdFormat<nind_remain>>::type
    minor_dims() const {
        constexpr auto nind_diff = nind-nind_remain;
        std::array<size_t, nind_remain> shape;
        std::copy(m_shape.cbegin()+nind_diff, m_shape.cend(), shape.begin());
        std::array<std::string, nind-nind_diff> dim_names;
        std::copy(m_dim_names.cbegin()+nind_diff, m_dim_names.cend(), dim_names.begin());
        return {shape, dim_names};
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

    /**
     * @return
     *  the shape array converted to a vector instance
     */
    defs::inds shape_vector() const {
        defs::inds tmp;
        tmp.assign(m_shape.cbegin(), m_shape.cend());
        return tmp;
    }

    const std::array<std::string, nind>& dim_names() const {
        return m_dim_names;
    }
    /**
     * @return
     *  the dim_names array converted to a vector instance
     */
    std::vector<std::string> dim_names_vector() const {
        std::vector<std::string> tmp;
        tmp.assign(m_dim_names.cbegin(), m_dim_names.cend());
        return tmp;
    }

    /**
     * NdFormat maps multidimensional indices to a single, contiguous integer index
     * @param inds
     *  array of multidimensional indices
     * @return
     *  flat index
     */
    size_t flatten(std::array<size_t, nind> inds) const {
        size_t i = 0ul;
        for (size_t iind=0ul; iind<nind; ++iind) i+=inds[iind]*m_strides[iind];
        return i;
    }

    /*
     * in addition to flattening a single array specifying the element completely, it is possible to flatten
     * multidimensional indices hierarchically. this is useful when the NdFormat is interpreted as the composition
     * of multiple multidimensional structures.
     *
     * suppose we have a 3-dimensional format which actually represents a 1-dimensional array of matrix indices. we may
     * want in this case to identify elements first by the matrix containing them, then by the pair of indices within
     * the identified matrix.
     *
     * This amounts to splitting up the inner product defining the flattening map into two parts
     */
    // TODO hierarchical flattening

    template<size_t nmajor>
    size_t combine(const size_t& iflat_major, const size_t& iflat_minor){
        static_assert(nmajor>0, "major flat index must correspond to non-zero number of dimensions");
        ASSERT(iflat_major < major_dims<nmajor>().nelement());
        auto stride = m_strides[nmajor-1];
        ASSERT(iflat_minor < stride);
        return iflat_major*stride+iflat_minor;
    }

    template<size_t nminor>
    size_t flatten(std::array<size_t, nind - nminor> major, std::array<size_t, nminor> minor) const {
        size_t iflat = 0ul;
        for (size_t i = 0ul; i < nind - nminor; ++i) {
            ASSERT(major[i]<m_shape[i]);
            iflat+= major[i] * m_strides[i];
        }
        for (size_t i = 0ul; i < nminor; ++i) {
            const auto j = i + nind - nminor;
            ASSERT(minor[i]<m_shape[j]);
            iflat+= minor[i] * m_strides[j];
        }
        return iflat;
    }

    template<size_t nminor>
    typename std::enable_if<nminor!=0, size_t>::type
    flatten(const size_t& iflat_major, const size_t& iflat_minor) const {
        return iflat_major*m_strides[nind-nminor-1]+iflat_minor;
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

    /**
     * helper function to store the strides assuming the m_shape member has already been set
     */
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
struct NdEnumeration : NdFormat<nind>{
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

    NdEnumeration(const NdFormat<nind>& format): NdFormat<nind>(format), m_inds(make_inds(format)){}

    const std::array<size_t, nind>& operator[](const size_t& i){
        return m_inds[i];
    }

};


#endif //M7_NDFORMAT_H