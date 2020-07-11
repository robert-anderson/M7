//
// Created by rja on 05/07/2020.
//

#ifndef M7_DISTRIBUTEDACCUMULATION_H
#define M7_DISTRIBUTEDACCUMULATION_H

#include "Distributed.h"

template<typename T, typename delta_T=T>
class DistributedAccumulation : public Distributed<T> {
public:
    using Distributed<T>::local;
    using Distributed<T>::mpi_sum;
    using Distributed<T>::operator=;
    using Distributed<T>::operator+=;
    Distributed<delta_T> m_delta;

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
