//
// Created by Robert John Anderson on 2020-02-13.
//

#ifndef M7_MUTEXVECTOR_H
#define M7_MUTEXVECTOR_H

#include <vector>
#include <omp.h>

class Mutex {
    omp_lock_t& m_lock;
    const size_t m_index;
public:
    Mutex(omp_lock_t &lock, const size_t &index);
    virtual ~Mutex();
    size_t index() const;
};


class MutexVector {
    std::vector<omp_lock_t> m_mutex;
public:
    MutexVector(size_t n);
    ~MutexVector();
    size_t size() const;
    void resize(const size_t &n);
    void grow(const size_t &n_add);
    Mutex get(const size_t &i);
};


#endif //M7_MUTEXVECTOR_H
