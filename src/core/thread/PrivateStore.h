//
// Created by Robert John Anderson on 2020-04-12.
//

#ifndef M7_PRIVATESTORE_H
#define M7_PRIVATESTORE_H

#include "AlignedAllocator.h"
#include "omp.h"

template<typename T>
class PrivateStore {
    struct alignas(defs::cache_line_size) aligned_T {
        T v;
    };
    static_assert(alignof(aligned_T)%defs::cache_line_size==0, "non-cache-aligned type");
    const size_t m_nthread;
    std::vector<aligned_T, AlignedAllocator<aligned_T, defs::cache_line_size>> m_data;

public:
    template<typename ...Args>
    PrivateStore(Args&&... construct_args):
    m_nthread(omp_get_max_threads()){
        for (size_t i = 0ul; i<m_nthread; ++i){
            m_data.push_back(aligned_T{T(std::forward<Args>(construct_args)...)});
        }
        // check alignment of first element
        ASSERT(((size_t)((void*)(m_data.data())))%defs::cache_line_size==0)
        // check that adjacent elements are properly spaced
        ASSERT(((size_t)((void*)(m_data.data()))-(size_t)((void*)(m_data.data()+1)) )%defs::cache_line_size==0)
        // check that adjacent elements are sufficiently spaced to contain type
        ASSERT(((size_t)((void*)(m_data.data()))-(size_t)((void*)(m_data.data()+1)) )>=sizeof(T))
    }

    PrivateStore():PrivateStore(1){}

    T& get(){
        auto const ithread = omp_get_thread_num();
        return m_data[ithread].v;
    }

    T* get_ptr(){
        auto const ithread = omp_get_thread_num();
        return (m_data.data()+ithread)->v;
    }

    const size_t& nthread() const {return m_nthread;}
};


#endif //M7_PRIVATESTORE_H
