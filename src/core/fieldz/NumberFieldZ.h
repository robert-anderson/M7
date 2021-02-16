//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELDZ_H
#define M7_NUMBERFIELDZ_H

#include "NdFieldBaseZ.h"

template<typename T, size_t nind_item, size_t nind_element>
struct NumberFieldBaseZ : FullyFormattedFieldBaseZ<T, nind_item, nind_element> {
    using FieldBaseZ::m_item_size;
    using FieldBaseZ::m_size;
    using FieldBaseZ::begin;
    using FieldBaseZ::end;
    typedef FullyFormattedFieldBaseZ<T, nind_item, nind_element> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    using base_t::nelement_all;

    NumberFieldBaseZ(std::array<size_t, nind_element> shape) :
            base_t(sizeof(T)*NdFormat<nind_element>(shape).nelement(), shape){}

    NumberFieldBaseZ &operator=(const T &v) {
        std::fill((T*)begin(), (T*)end(), v);
        return *this;
    }

    NumberFieldBaseZ &operator=(const std::vector<T> &v) {
        ASSERT(v.size() == nelement_all());
        std::memcpy(begin(), v.data(), m_size);
        return *this;
    }

    T& operator[](const std::pair<size_t, size_t>& pair){
        const auto& iitem = pair.first;
        const auto& ielement = pair.second;
        return ((T *) begin(iitem))[ielement];
    }

    const T& operator[](const std::pair<size_t, size_t>& pair) const {
        const auto& iitem = pair.first;
        const auto& ielement = pair.second;
        return ((T *) begin(iitem))[ielement];
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string tmp;
        if (nind_element) tmp += "[";
        for (size_t ielement = 0ul; ielement<m_element_format.nelement(); ++ielement)
            tmp+=std::to_string((*this)[{iitem, ielement}]) + " ";
        if (nind_element) tmp += "]";
        return tmp;
    }
};

/*
 * The following definitions are to make accesses prettier.
 */


template<typename T, size_t nind_item, size_t nind_element>
struct NumberFieldZ : NumberFieldBaseZ<T, nind_item, nind_element> {
    typedef NumberFieldBaseZ<T, nind_item, nind_element> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    NumberFieldZ(std::array<size_t, nind_element> shape) : NumberFieldBaseZ<T, nind_item, nind_element>(shape){}

    T &operator()(const std::array<size_t, nind_item>& inds, const std::array<size_t, nind_element>& einds) {
        return (*this)[{m_item_format->flatten(inds), m_element_format.flatten(einds)}];
    }

    const T &operator()(const std::array<size_t, nind_item>& inds, const std::array<size_t, nind_element>& einds) const {
        return (*this)[{m_item_format->flatten(inds), m_element_format.flatten(einds)}];
    }
};


template<typename T, size_t nind_item>
struct NumberFieldZ<T, nind_item, 0ul> : NumberFieldBaseZ<T, nind_item, 0ul> {
    typedef NumberFieldBaseZ<T, nind_item, 0ul> base_t;
    using base_t::m_item_format;
    NumberFieldZ() : NumberFieldBaseZ<T, nind_item, 0ul>({}){}

    T &operator()(const std::array<size_t, nind_item>& inds) {
        return (T *) begin(m_item_format->flatten(inds));
    }

    const T &operator()(const std::array<size_t, nind_item>& inds) const {
        return (const T *) begin(m_item_format->flatten(inds));
    }
};

template<typename T, size_t nind_element>
struct NumberFieldZ<T, 0ul, nind_element> : NumberFieldBaseZ<T, 0ul, nind_element> {
    typedef NumberFieldBaseZ<T, 0ul, nind_element> base_t;
    using base_t::m_element_format;
    using base_t::begin;
    NumberFieldZ(std::array<size_t, nind_element> shape) : NumberFieldBaseZ<T, 0ul, nind_element>(shape){}

    T &operator()(const std::array<size_t, nind_element>& einds) {
        return ((T *) begin())[m_element_format.flatten(einds)];
    }

    const T &operator()(const std::array<size_t, nind_element>& einds) const {
        return ((const T *) begin())[m_element_format.flatten(einds)];
    }
};

template<typename T>
struct NumberFieldZ<T, 0ul, 0ul> : NumberFieldBaseZ<T, 0ul, 0ul> {
    typedef NumberFieldBaseZ<T, 0ul, 0ul> base_t;
    using base_t::m_element_format;
    using base_t::begin;
    NumberFieldZ() : NumberFieldBaseZ<T, 0ul, 0ul>({}){}

    T &operator()() {
        return (T *) begin();
    }

    const T &operator()() const {
        return (const T *) begin();
    }

};

template <size_t nind_item>
using BosonOnvFieldZ = NumberFieldZ<uint8_t, nind_item, 1ul>;


#endif //M7_NUMBERFIELDZ_H
