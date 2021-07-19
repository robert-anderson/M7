//
// Created by rja on 12/10/2020.
//

#ifndef M7_NDACCESSOR_H
#define M7_NDACCESSOR_H

#include <cstring>
#include <functional>
#include "NdFormat.h"
#include "src/core/parallel/MPIWrapper.h"

template <typename T, size_t nind>
struct NdAccessor {
    T* m_data;
    const NdFormat<nind>& m_format;

    template<typename ...Args>
    NdAccessor(T* data, const NdFormat<nind>& format):
    m_data{data}, m_format(format){}

    NdAccessor(const NdAccessor &other): m_data(other.m_data), m_format(other.m_format){}

    NdAccessor& operator=(const NdAccessor &other){
        if (&other!=this){
            ASSERT(other.m_format==m_format);
            ASSERT(other.nelement()==nelement());
            std::memcpy(m_data, other.m_data, nelement()*sizeof(T));
        }
        return *this;
    }

    template<typename U>
    NdAccessor& operator=(const std::vector<U> &v){
        ASSERT(v.size()==nelement());
        for (size_t i = 0ul; i < nelement(); ++i) (*this)(i) = v[i];
        return *this;
    }

    NdAccessor &operator=(const std::vector<T> &v) {
        ASSERT(v.size() == nelement());
        std::memcpy(
                reinterpret_cast<void*>(m_data),
                reinterpret_cast<void*>(v.data()),
                nelement()*sizeof(T));
        return *this;
    }

    T& operator[](const size_t& ielement){
        ASSERT(ielement<m_format.nelement());
        return m_data[ielement];
    }

    const T& operator[](const std::array<size_t, nind>& inds) const {
        return m_data[m_format.flatten(inds)];
    }

    T& operator[](const std::array<size_t, nind>& inds) {
        return m_data[m_format.flatten(inds)];
    }

    template<typename ...Args>
    T& operator()(Args... inds){
        return m_data[m_format.flatten(inds...)];
    }

    template<typename ...Args>
    const T& operator()(Args... inds) const {
        return m_data[m_format.flatten(inds...)];
    }

    void bcast(size_t iroot=0ul) {
        mpi::bcast(m_data, nelement(), iroot);
    }

    size_t nelement() const {
        return m_format.nelement();
    }

    void clear() {
        std::memset(reinterpret_cast<void*>(m_data), 0, sizeof(T)*nelement());
    }
};

/*
 * An NdAccessor that owns its own NdFormat, rather than referencing it from outside
 */
template <typename T, size_t nind>
struct FormattedNdAccessor : private NdFormat<nind>, public NdAccessor<T, nind>{
    template<typename ...Args>
    FormattedNdAccessor(T* data, NdFormat<nind> format):
    NdFormat<nind>(format.shape()),
    NdAccessor<T, nind>(data, *this){}
};

#endif //M7_NDACCESSOR_H
