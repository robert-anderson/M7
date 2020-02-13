//
// Created by Robert John Anderson on 2020-02-13.
//

#ifndef M7_MUTEXVECTOR_H
#define M7_MUTEXVECTOR_H

#include <vector>
#include <omp.h>


class MutexVector {
    std::vector<omp_lock_t> m_mutex;
public:
    MutexVector(size_t n);
    ~MutexVector();
    void grow(const size_t &n_add);
    void acquire_lock(size_t i);
    void release_lock(size_t i);
};


#endif //M7_MUTEXVECTOR_H
