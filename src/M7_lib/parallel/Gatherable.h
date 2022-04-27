//
// Created by Robert J. Anderson on 12/07/2020.
//

#ifndef M7_GATHERABLE_H
#define M7_GATHERABLE_H

#include <vector>
#include "Distributed.h"

template<typename T>
class Gatherable : public Distributed<T> {
protected:
    std::vector<T> m_gathered{}; // result of MPI gather
public:
    using Distributed<T>::m_local;
    using Distributed<T>::operator=;
    Gatherable() : m_gathered(mpi::nrank()) {}

    const std::vector<T> &gathered() {
        return m_gathered;
    }

    const std::vector<T> &mpi_gather() {
        mpi::all_gather(m_local, m_gathered);
        return m_gathered;
    }

};


#endif //M7_GATHERABLE_H
