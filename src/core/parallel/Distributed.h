//
// Created by rja on 15/06/2020.
//

#ifndef M7_DISTRIBUTED_H
#define M7_DISTRIBUTED_H

#include "MPIWrapper.h"
#include "src/core/thread/Atomic.h"

template<typename T>
class Distributed {
    static_assert(std::is_arithmetic<T>::value, "template type must be arithmetic.");
protected:
    T m_local = 0; // value on this process
    T m_reduced = std::numeric_limits<T>::max(); // cached result of last all_reduce operation

public:

    Distributed(){
        if (omp_get_level()) throw std::runtime_error("Distributed variable must be declared outside parallel region");
    }

    virtual Distributed<T>& operator=(const T& rhs){
        m_local = rhs;
        m_reduced = std::numeric_limits<T>::max();
        return *this;
    }

    virtual Distributed<T>& operator+=(const T& rhs){
        m_local += rhs;
        m_reduced = std::numeric_limits<T>::max();
        return *this;
    }

    T& local(){
        return m_local;
    }

    T& reduced(){
#ifndef DNDEBUG
        // check that a reduction method was applied
        //ASSERT(m_reduced!=std::numeric_limits<T>::max());
#endif
        return m_reduced;
    }

    T& mpi_sum(){
        m_reduced = mpi::all_sum(m_local);
        return m_reduced;
    }

    T& mpi_max(){
        m_reduced = mpi::all_max(m_local);
        return m_reduced;
    }

    T& mpi_min(){
        m_reduced = mpi::all_min(m_local);
        return m_reduced;
    }

    T& mpi_bcast(size_t irank){
        m_reduced = m_local;
        mpi::bcast(m_reduced, irank);
    }
};


#endif //M7_DISTRIBUTED_H
