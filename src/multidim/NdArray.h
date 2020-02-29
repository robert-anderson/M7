//
// Created by Robert John Anderson on 2020-02-22.
//

#ifndef M7_NDARRAY_H
#define M7_NDARRAY_H

#include <cstddef>
#include "Indexer.h"

/*
 * NdArray is a fully and partially subscriptable multidimensional array
 * implementation. An instance may manage its own memory, or constitute a
 * view on another instance, typically as a sub-array, e.g. a particular
 * row of a matrix (NdArray<2>).
 */

template<typename T, size_t nind>
class NdArray {
protected:
    const Indexer<nind> m_indexer;
    std::vector<T> m_data_internal;
    T* m_data = nullptr;

public:
    template<typename ...Args>
    NdArray(const size_t &first, Args ...shape) : m_indexer(first, shape...),
                             m_data_internal(nelement(), 0), m_data(m_data_internal.data()) {}

    template<typename ...Args>
    NdArray(T* data, const std::array<size_t, nind> &shape) : m_indexer(shape), m_data(data){}

    template<typename ...Args>
    T *view(Args ...inds) const {
        return m_data+m_indexer.get(inds...);
    }

    template<typename ...Args>
    NdArray<T, nind-sizeof...(Args)> subarray(Args ...subs){
        std::array<size_t, nind-sizeof...(subs)> shape;
        std::copy(m_indexer.shape().begin()+sizeof...(subs),
            m_indexer.shape().end(), shape.begin());
        return NdArray<T, nind-sizeof...(subs)>(m_data+m_indexer.get_sub(subs...), shape);
    }

    void operator=(NdArray &src) {
        assert(m_indexer.shape() == src.m_shape);
        memcpy(m_data, src.m_data, nelement() * sizeof(T));
    }

    void operator=(std::vector<T> &src) {
        assert(nelement() == src.size());
        memcpy(m_data, src.data(), nelement() * sizeof(T));
    }

    void operator*=(const T& factor){
        for (auto i=m_data; i!=m_data+nelement(); ++i) *i*=factor;
    }

    void operator/=(const T& factor){
        for (auto i=m_data; i!=m_data+nelement(); ++i) *i/=factor;
    }

    const std::array<size_t, nind> &shape() const {
        return m_indexer.shape();
    }

    const std::array<size_t, nind> &strides() const {
        return m_indexer.strides();
    }

    size_t nelement() const {
        return m_indexer.nelement();
    }


};


#endif //M7_NDARRAY_H
