//
// Created by rja on 23/06/2020.
//

#ifndef M7_HYBRID_H
#define M7_HYBRID_H

#include "src/core/parallel/Distributed.h"
#include "src/core/thread/PrivateStore.h"

template<typename T>
class Hybrid : public Distributed<T> {
    PrivateStore<T> m_thread;
    using Distributed<T>::m_local;
public:
    T& thread(){
        ASSERT(omp_get_level())
        return m_thread.get();
    }

    void zero_thread_values(){
        ASSERT(!omp_get_level())
        return m_thread.zero();
    }

    Hybrid<T>& operator=(const T& rhs) override {
        if (omp_get_level()) thread() = rhs;
        else Distributed<T>::operator=(rhs);
        return *this;
    }

    Hybrid<T>& operator+=(const T& rhs) override {
        if (omp_get_level()) thread() += rhs;
        else Distributed<T>::operator+=(rhs);
        return *this;
    }

    T thread_sum(){
        ASSERT(!omp_get_level())
        return m_thread.reduce_sum();
    }

    T& put_thread_sum(){
        ASSERT(!omp_get_level())
        m_local=thread_sum();
        zero_thread_values();
        return m_local;
    }

    T& add_thread_sum(){
        ASSERT(!omp_get_level())
        m_local+=thread_sum();
        zero_thread_values();
        return m_local;
    }

    using Distributed<T>::mpi_sum;
    T& sum(){
        m_local = thread_sum();
        return mpi_sum();
    }

    T thread_max(){
        ASSERT(!omp_get_level())
        return m_thread.reduce_max();
    }

    using Distributed<T>::mpi_max;
    T& max(){
        m_local = thread_max();
        return mpi_max();
    }

    T& thread_min(){
        ASSERT(!omp_get_level())
        return m_thread.reduce_min();
    }

    using Distributed<T>::mpi_min;
    T& min(){
        m_local = thread_max();
        return mpi_min();
    }
};


#endif //M7_HYBRID_H
