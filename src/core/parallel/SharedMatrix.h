//
// Created by Robert John Anderson on 2020-08-10.
//

#ifndef M7_SHAREDMATRIX_H
#define M7_SHAREDMATRIX_H

#include "SharedArray.h"

template<typename T>
class SharedMatrix : public SharedArray<T> {
    const size_t m_nrow, m_ncol;

public:
    SharedMatrix(size_t nrow, size_t ncol):
        SharedArray<T>(nrow*ncol), m_nrow(nrow), m_ncol(ncol){}

    void set(const size_t &irow, const size_t& icol, const T& v) {
        // element-modifying access should only take place on the root rank
        ASSERT(irow<m_nrow)
        ASSERT(icol<m_ncol)
        SharedArray<T>::set(irow*m_ncol+icol, v);
    }

    const T& get(const size_t &irow, const size_t& icol) const {
        ASSERT(irow<m_nrow)
        ASSERT(icol<m_ncol)
        return SharedArray<T>::get(irow*m_ncol+icol);
    }

};


#endif //M7_SHAREDMATRIX_H
