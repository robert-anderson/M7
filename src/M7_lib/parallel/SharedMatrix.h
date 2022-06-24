//
// Created by Robert John Anderson on 2020-08-10.
//

#ifndef M7_SHAREDMATRIX_H
#define M7_SHAREDMATRIX_H

#include "SharedArray.h"

template<typename T>
class SharedMatrix : public SharedArray<T> {
    const uint_t m_nrow, m_ncol;

public:
    SharedMatrix(uint_t nrow, uint_t ncol):
        SharedArray<T>(nrow*ncol), m_nrow(nrow), m_ncol(ncol){}

    void set(const uint_t &irow, const uint_t& icol, const T& v) {
        // element-modifying access should only take place on the root rank
        ASSERT(irow<m_nrow)
        ASSERT(icol<m_ncol)
        SharedArray<T>::set(irow*m_ncol+icol, v);
    }

    const T& get(const uint_t &irow, const uint_t& icol) const {
        ASSERT(irow<m_nrow)
        ASSERT(icol<m_ncol)
        return SharedArray<T>::operator[](irow*m_ncol+icol);
    }

};


#endif //M7_SHAREDMATRIX_H
