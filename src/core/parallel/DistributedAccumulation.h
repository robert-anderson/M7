//
// Created by rja on 05/07/2020.
//

#ifndef M7_DISTRIBUTEDACCUMULATION_H
#define M7_DISTRIBUTEDACCUMULATION_H

#include "Hybrid.h"

template<typename T, typename delta_T=T>
class DistributedAccumulation : public Distributed<T> {
public:
    using Distributed<T>::local;
    using Distributed<T>::mpi_sum;
    using Distributed<T>::operator=;
    using Distributed<T>::operator+=;
    Hybrid<delta_T> m_delta;

    delta_T& reduce_delta(){
        m_delta.put_thread_sum();
        return m_delta.mpi_sum();
    }

    T& accumulate(){
        local() += m_delta.local();
        return mpi_sum();
    }
};


#endif //M7_DISTRIBUTEDACCUMULATION_H
