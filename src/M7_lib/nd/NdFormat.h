//
// Created by rja on 02/10/2020.
//

#ifndef M7_NDFORMAT_H
#define M7_NDFORMAT_H

#include <M7_lib/enumerator/Enumerator.h>
#include <M7_lib/defs.h>
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
struct NdFormat : NdFormatBase {
    /**
     * the extents of each dimension
     */
    const std::array<size_t, nind> m_shape;
    /**
     * elements of the stride are computed as products of elements of the shape. this allows for the conversion of an
     * array containing multidimensional indices to a single integer ("flat" index) via a scalar product
     */
    const std::array<size_t, nind> m_strides;
    /**
     * label each dimension
     */
    const std::array<std::string, nind> m_dim_names;
    /**
     * the number of elements in the indexing scheme. the largest value a flat index can take in this scheme is one less
     * than m_nelement
     */
    const size_t m_nelement;

    /**
     * helper function to compute the strides assuming the m_shape member has already been set
     */
    std::array<size_t, nind> make_strides() const {
        std::array<size_t, nind> tmp{};
        if (!nind) return {};
        tmp.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            tmp[nind - i] = tmp[nind - i + 1] * m_shape[nind - i + 1];
        }
        return tmp;
    }

    /**
     * helper function to decide the number of elements assuming the m_shape and m_stride members have already been set
     */
    size_t make_nelement() const {
        return nind ? m_shape.front()*m_strides.front() : 1ul;
    }

    std::array<std::string, nind> make_default_dim_names() const {
        std::array<std::string, nind> dim_names;
        for (size_t i=0ul; i!=nind; ++i) dim_names[i] = "dim_"+std::to_string(i);
        return dim_names;
    }


public:
    NdFormat(): m_shape(), m_strides(), m_dim_names(), m_nelement(make_nelement()){
        static_assert(!nind, "This ctor is only valid in the scalar case");
    }

    NdFormat(const std::array<size_t, nind>& shape, const std::array<std::string, nind>& dim_names):
        m_shape(shape), m_strides(make_strides()), m_dim_names(dim_names), m_nelement(make_nelement()){
        ASSERT(m_nelement!=~0ul);
    }

    NdFormat(const std::array<size_t, nind>& shape): NdFormat(shape, make_default_dim_names()){}

    /**
     * construct a scheme where all extents are the same
     * @param extent
     *  value to copy to all elements of the shape
     */
    NdFormat(size_t extent): NdFormat(array_utils::filled<size_t, nind>(extent)){}

    NdFormat(const NdFormat<nind>& other) : NdFormat(other.m_shape, other.m_dim_names){}

    bool operator==(const NdFormat<nind> &other) const{
        for (size_t iind=0ul; iind!=nind; ++iind){
            if (other.m_shape[iind]!=m_shape[iind]) return false;
        }
        return true;
    }

    operator const std::array<size_t, nind>&() const{
        return m_shape;
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

    /**
     * @return
     *  the shape array converted to a vector instance
     */
    defs::inds shape_vector() const {
        return array_utils::to_vector(m_shape);
    }

    /**
     * @return
     *  the dim_names array converted to a vector instance
     */
    std::vector<std::string> dim_names_vector() const {
        return array_utils::to_vector(m_dim_names);
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
        for (size_t iind=0ul; iind!=nind; ++iind) {
            DEBUG_ASSERT_LT(inds[iind], m_shape[iind], "index OOB");
            i+=inds[iind]*m_strides[iind];
        }
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
        ASSERT(iflat_major < major_dims<nmajor>().m_nelement);
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

    template<typename T>
    void decode_flat(const size_t& iflat, std::array<T, nind>& inds) const {
        static_assert(std::is_integral<T>::value, "index type must be integral");
        size_t remainder = iflat;
        for (size_t i=0ul; i!=nind; ++i){
            auto& ind = inds[i];
            ind = remainder/m_strides[i];
            remainder-=ind*m_strides[i];
        }
    }

    std::string to_string() const override {
        if (!nind) return "scalar";
        std::string tmp;
        for (size_t i=0ul; i!=nind; ++i) {
            tmp+=m_dim_names[i]+" ("+std::to_string(m_shape[i])+") ";
        }
        return tmp;
    }


private:

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
        for (size_t i=0ul; i!=nind_spec; ++i) {
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
        std::vector<std::array<size_t, nind>> out(format.m_nelement);
        enums::PermutationsWithRepetition e(format.shape_vector());
        size_t iinds = 0ul;
        while (e.next()){
            for (size_t i=0; i!=nind; ++i) out[iinds][i] = e[i];
            ++iinds;
        }
        ASSERT(iinds==format.m_nelement);
        return out;
    }

    NdEnumeration(const NdFormat<nind>& format): NdFormat<nind>(format), m_inds(make_inds(format)){}

    const std::array<size_t, nind>& operator[](const size_t& i){
        return m_inds[i];
    }

};


#endif //M7_NDFORMAT_H