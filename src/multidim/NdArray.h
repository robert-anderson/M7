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
    std::vector<T> m_data;

public:
    template<typename ...Args>
    NdArray(const size_t &first, Args ...shape) : m_indexer(first, shape...),
                             m_data(nelement(), 0) {}


    template<typename ...Args>
    NdArray(const std::array<size_t, nind> &shape) : m_indexer(shape){}

    template<typename ...Args>
    T *view(Args ...inds) const {
        return (T*)m_data.data()+m_indexer.get(inds...);
    }

    template<typename ...Args>
    auto subarray(const size_t &first, Args ...subs){
        std::array<size_t, nind-sizeof...(subs)-1> shape;
        std::copy(m_indexer.shape().begin()+sizeof...(subs)-1,
            m_indexer.shape().end(), shape.begin());
        return NdArray<T, nind-sizeof...(subs)-1>(shape);
    }

    void operator=(NdArray &src) {
        assert(Indexer<nind>::m_shape == src.m_shape);
        memcpy(m_data.data(), src.m_data, nelement() * sizeof(T));
    }

    void operator=(std::vector<T> &src) {
        assert(Indexer<nind>::nelement() == src.size());
        memcpy(m_data.data(), src.data(), nelement() * sizeof(T));
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
