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

    void zero() {
        std::memset((void*)m_data, 0, sizeof(T)*nelement());
    }
};

template <typename T, size_t nind>
struct FormattedNdAccessor : private NdFormat<nind>, public NdAccessor<T, nind>{
    template<typename ...Args>
    FormattedNdAccessor(T* data, NdFormat<nind> format):
    NdFormat<nind>(format.shape()),
    NdAccessor<T, nind>(data, *this){}
};

#endif //M7_NDACCESSOR_H
