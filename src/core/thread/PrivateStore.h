//
// Created by Robert John Anderson on 2020-04-12.
//

#ifndef M7_PRIVATESTORE_H
#define M7_PRIVATESTORE_H

#include "AlignedAllocator.h"
#include "omp.h"

template<typename T>
class PrivateStore {
    //static_assert(alignof(T)%defs::cache_line_size==0, "Cannot store non-cache-aligned type in PrivateStore");
    const size_t m_nthread;
    const size_t m_nelement_per_thread;
    std::vector<T, AlignedAllocator<T, defs::cache_line_size>> m_data;

public:
    PrivateStore(size_t nelement_per_thread, const T& t):
    m_nthread(omp_get_max_threads()), m_nelement_per_thread(nelement_per_thread),
    m_data(m_nthread*nelement_per_thread, t){}

    PrivateStore():PrivateStore(1, 0){}

    T& get(const size_t &ielement=0){
        auto const ithread = omp_get_thread_num();
        return m_data[m_nelement_per_thread*ithread+ielement];
    }

    const size_t& nthread() const {return m_nthread;}
};


#endif //M7_PRIVATESTORE_H
