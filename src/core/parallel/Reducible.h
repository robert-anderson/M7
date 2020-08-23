//
// Created by rja on 15/06/2020.
//

#ifndef M7_REDUCIBLE_H
#define M7_REDUCIBLE_H

#include "MPIWrapper.h"
#include "Distributed.h"

template<typename T>
class Reducible : public Distributed<T> {
protected:
    std::pair<T, size_t> m_reduced {std::numeric_limits<T>::max(),
                                    std::numeric_limits<size_t>::max()}; // cached result of last all_reduce operation
public:
    using Distributed<T>::m_local;
    using Distributed<T>::operator=;

    Reducible(){}

    const T& reduced() const {
        return m_reduced.first;
    }

    const size_t& irank() const {
        return m_reduced.second;
    }

    const T& mpi_sum(){
        m_reduced.first = mpi::all_sum(m_local);
        return m_reduced.first;
    }

    const T& mpi_max(){
        m_reduced.first = mpi::all_max(m_local);
        return m_reduced.first;
    }

    const T& mpi_min(){
        m_reduced.first = mpi::all_min(m_local);
        return m_reduced.first;
    }

    const T& mpi_bcast(size_t irank){
        m_reduced.first = m_local;
        mpi::bcast(m_reduced.first, irank);
        return m_reduced.first;
    }

    const size_t& mpi_minloc(){
        mpi::all_minloc(m_local, m_reduced);
        return m_reduced.second;
    }

    const size_t& mpi_maxloc(){
        mpi::all_maxloc(m_local, m_reduced);
        return m_reduced.second;
    }

};


#endif //M7_REDUCIBLE_H
