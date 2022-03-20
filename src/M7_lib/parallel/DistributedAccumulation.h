//
// Created by rja on 05/07/2020.
//

#ifndef M7_DISTRIBUTEDACCUMULATION_H
#define M7_DISTRIBUTEDACCUMULATION_H

#include "Reducible.h"

template<typename T, typename delta_T=T>
class DistributedAccumulation : public Reducible<T> {
public:
    using Reducible<T>::local;
    using Reducible<T>::mpi_sum;
    using Reducible<T>::operator=;
    using Reducible<T>::operator+=;
    Reducible<delta_T> m_delta;

    const delta_T& reduce_delta(){
        return m_delta.mpi_sum();
    }

    const T& accumulate(){
        local() += m_delta.local();
        m_delta.reset();
        return mpi_sum();
    }
};


#endif //M7_DISTRIBUTEDACCUMULATION_H
