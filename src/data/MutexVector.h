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
    size_t size() const;
    void resize(const size_t &n);
    void grow(const size_t &n_add);
    void acquire_lock(const size_t &i);
    void release_lock(const size_t &i);
};


#endif //M7_MUTEXVECTOR_H
