//
// Created by Robert John Anderson on 2020-02-13.
//

#include <assert.h>
#include <iostream>
#include "MutexVector.h"

MutexVector::MutexVector(size_t n) {
    resize(n);
}

void MutexVector::resize(const size_t &n) {
    for (size_t i = 0; i < m_mutex.size(); ++i) omp_destroy_lock(m_mutex.data() + i);
    m_mutex.resize(n);
    for (size_t i = 0; i < m_mutex.size(); ++i) omp_init_lock(m_mutex.data() + i);
}

void MutexVector::grow(const size_t &n_add) {
    resize(m_mutex.size() + n_add);
}

void MutexVector::acquire_lock(const size_t &i) {
    assert(i<size());
    omp_set_lock(m_mutex.data() + i);
}

void MutexVector::release_lock(const size_t &i) {
    assert(i<size());
    omp_unset_lock(m_mutex.data() + i);
}

MutexVector::~MutexVector() {
    for (auto i{0ul}; i < m_mutex.size(); ++i) omp_destroy_lock(m_mutex.data() + i);
}

size_t MutexVector::size() const {
    return m_mutex.size();
}
