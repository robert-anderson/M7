//
// Created by Robert John Anderson on 2020-02-13.
//

#include "MutexVector.h"

MutexVector::MutexVector(size_t n) :m_mutex(n){
    for (auto i{0ul}; i<n; ++i) omp_init_lock(m_mutex.data()+i);
}

void MutexVector::grow(const size_t &n_add) {
    m_mutex.resize(m_mutex.size()+n_add);
    for (size_t i = m_mutex.size()-n_add; i<m_mutex.size(); ++i){
        omp_init_lock(m_mutex.data()+i);
    }
}

void MutexVector::acquire_lock(size_t i) {
    omp_set_lock(m_mutex.data()+i);
}

void MutexVector::release_lock(size_t i) {
    omp_unset_lock(m_mutex.data()+i);
}

MutexVector::~MutexVector() {
    for (auto i{0ul}; i<m_mutex.size(); ++i) omp_destroy_lock(m_mutex.data()+i);
}
