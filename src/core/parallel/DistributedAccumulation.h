//
// Created by rja on 05/07/2020.
//

#ifndef M7_DISTRIBUTEDACCUMULATION_H
#define M7_DISTRIBUTEDACCUMULATION_H

#include "Reducable.h"

template<typename T, typename delta_T=T>
class DistributedAccumulation : public Reducable<T> {
public:
    using Reducable<T>::local;
    using Reducable<T>::mpi_sum;
    using Reducable<T>::operator=;
    using Reducable<T>::operator+=;
    Reducable<delta_T> m_delta;

    delta_T& reduce_delta(){
        return m_delta.mpi_sum();
    }

    T& accumulate(){
        local() += m_delta.local();
        m_delta.reset();
        return mpi_sum();
    }
};


#endif //M7_DISTRIBUTEDACCUMULATION_H
