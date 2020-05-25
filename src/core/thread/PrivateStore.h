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
    std::vector<T, AlignedAllocator<T, defs::cache_line_size>> m_data;

public:
    template<typename ...Args>
    PrivateStore(Args... construct_args):
    m_nthread(omp_get_max_threads())
    {
        for (size_t i = 0ul; i<m_nthread; ++i){
            m_data.emplace_back(construct_args...);
        }
    }

    PrivateStore():PrivateStore(1){}

    T& get(){
        auto const ithread = omp_get_thread_num();
        return m_data[ithread];
    }

    const size_t& nthread() const {return m_nthread;}
};


#endif //M7_PRIVATESTORE_H
