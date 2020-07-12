//
// Created by rja on 15/06/2020.
//

#ifndef M7_REDUCABLE_H
#define M7_REDUCABLE_H

#include "MPIWrapper.h"
#include "Distributed.h"

template<typename T>
class Reducable : public Distributed<T> {
protected:
    T m_reduced = std::numeric_limits<T>::max(); // cached result of last all_reduce operation
public:
    using Distributed<T>::m_local;
    using Distributed<T>::operator=;

    Reducable(){}

    const T& reduced() const {
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
        return m_reduced;
    }
};


#endif //M7_REDUCABLE_H
