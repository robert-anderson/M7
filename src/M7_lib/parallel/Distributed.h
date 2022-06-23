//
// Created by Robert J. Anderson on 12/07/2020.
//

#ifndef M7_DISTRIBUTED_H
#define M7_DISTRIBUTED_H

#include "MPIWrapper.h"

template<typename T>
class Distributed {
    static_assert(datatype::is_arithmetic<T>(), "template type must be arithmetic.");
protected:
    T m_local{}; // value on this process

public:

    Distributed(){}

    operator T&() { return m_local; }
    operator T() const { return m_local; }

    Distributed<T>& operator=(const T& rhs){
        m_local = rhs;
        return *this;
    }

    void reset(){
        *this=0;
    }

    Distributed<T>& operator+=(const T& rhs){
        m_local += rhs;
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

};


#endif //M7_DISTRIBUTED_H
