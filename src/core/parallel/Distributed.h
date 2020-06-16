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
    T m_local = 0; // value on this process
    T m_reduced = std::numeric_limits<T>::max(); // cached result of last all_reduce operation

public:
    Distributed<T>& operator=(const T& rhs){
        m_local = rhs;
        m_reduced = std::numeric_limits<T>::max();
        return *this;
    }

    //operator T&() { return local; }

    T& local(){
        return m_local;
    }

    T& reduced(){
#ifndef DNDEBUG
        // check that a reduction method was applied
        ASSERT(m_reduced!=std::numeric_limits<T>::max());
#endif
        return m_reduced;
    }

    T& sum(){
        m_reduced = mpi::all_sum(m_local);
        return m_reduced;
    }

    T& max(){
        m_reduced = mpi::all_max(m_local);
        return m_reduced;
    }

    T& min(){
        m_reduced = mpi::all_min(m_local);
        return m_reduced;
    }

    T& bcast(size_t irank){
        m_reduced = m_local;
        mpi::bcast(m_reduced, irank);
    }

    Atomic<T> as_atomic() {
        return Atomic<T>(m_local);
    }
};


#endif //M7_DISTRIBUTED_H
