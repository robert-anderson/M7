//
// Created by rja on 15/06/2020.
//

#ifndef M7_DISTRIBUTED_H
#define M7_DISTRIBUTED_H

#include "MPIWrapper.h"

template<typename T>
class Distributed {
    static_assert(std::is_arithmetic<T>::value, "template type must be arithmetic.");
protected:
    T m_local = 0; // value on this process
    T m_reduced = std::numeric_limits<T>::max(); // cached result of last all_reduce operation

public:

    Distributed(){
    }

    operator T&() { return m_local; }
    //operator T() const { return m_local; }

    Distributed<T>& operator=(const T& rhs){
        m_local = rhs;
        m_reduced = std::numeric_limits<T>::max();
        return *this;
    }

    void reset(){
        *this=0;
    }

    Distributed<T>& operator+=(const T& rhs){
        m_local += rhs;
        m_reduced = std::numeric_limits<T>::max();
        return *this;
    }

    /*
    bool operator==(const T &rhs) const {
        return m_local == rhs;
    }

    bool operator!=(const T &rhs) const {
        return !(rhs == *this);
    }

    const Distributed<T> operator++(int) {
        m_local++;
        return *this;
    }

    const Distributed<T> operator--(int) {
        m_local--;
        return *this;
    }
     */

    T& local() {
        return m_local;
    }

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


#endif //M7_DISTRIBUTED_H
