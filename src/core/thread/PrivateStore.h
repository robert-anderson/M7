//
// Created by Robert John Anderson on 2020-04-12.
//

#ifndef M7_PRIVATESTORE_H
#define M7_PRIVATESTORE_H

#include "AlignedAllocator.h"
#include "omp.h"
#include <numeric>

template<typename T>
class PrivateStore {
    struct  aligned_T {
        T v;
    };
    //static_assert(alignof(aligned_T)%defs::cache_line_size==0, "non-cache-aligned type");
    const size_t m_nthread;
    std::vector<aligned_T, AlignedAllocator<aligned_T, defs::cache_line_size>> m_data;

    template< bool cond, typename U >
    using enable_if_t  = typename std::enable_if< cond, U >::type;

public:
    template<typename ...Args>
    PrivateStore(Args&&... construct_args):
    m_nthread(omp_get_max_threads()){
        for (size_t i = 0ul; i<m_nthread; ++i){
            m_data.push_back(aligned_T{T(std::forward<Args>(construct_args)...)});
        }
        return;
        // check alignment of first element
        ASSERT(((size_t)((void*)(m_data.data())))%defs::cache_line_size==0)
        // check that adjacent elements are properly spaced
        ASSERT(((size_t)((void*)(m_data.data()))-(size_t)((void*)(m_data.data()+1)) )%defs::cache_line_size==0)
        // check that adjacent elements are sufficiently spaced to contain type
        ASSERT(((size_t)((void*)(m_data.data()))-(size_t)((void*)(m_data.data()+1)) )>=sizeof(T))
    }

    T& get(){
        auto const ithread = omp_get_thread_num();
        ASSERT(ithread<m_nthread);
        return m_data[ithread].v;
    }

    const size_t& nthread() const {return m_nthread;}

    template<typename U = T>
    enable_if_t<true, U> reduce_land(){
        T val = 1;
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) val&= m_data[ithread].v;
        return val;
    }

    template<typename U = T>
    enable_if_t<true, U> reduce_lor(){
        T val = 0;
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) val|= m_data[ithread].v;
        return val;
    }

    template<typename U = T>
    enable_if_t<true, U> reduce_sum(){
        T val = 0;
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) val+= m_data[ithread].v;
        return val;
    }

    template<typename U = T>
    enable_if_t<true, U> reduce_prod(){
        T val = 1;
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) val*= m_data[ithread].v;
        return val;
    }

    template<typename U = T>
    enable_if_t<std::is_integral<U>::value, U> reduce_max(){
        T val = std::numeric_limits<T>::min();
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) {
            if (m_data[ithread].v>val) val = m_data[ithread].v;
        }
        return val;
    }

    template<typename U = T>
    enable_if_t<!std::is_integral<U>::value, U> reduce_max(){
        T val = std::numeric_limits<T>::min();
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) {
            auto tmp = std::abs(m_data[ithread].v);
            if (tmp>val) val = tmp;
        }
        return val;
    }

    template<typename U = T>
    enable_if_t<std::is_integral<U>::value, U> reduce_min(){
        T val = std::numeric_limits<T>::max();
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) {
            if (m_data[ithread].v<val) val = m_data[ithread].v;
        }
        return val;
    }

    template<typename U = T>
    enable_if_t<!std::is_integral<U>::value, U> reduce_min(){
        T val = std::numeric_limits<T>::min();
        for (size_t ithread=0ul; ithread<nthread(); ++ithread) {
            auto tmp = std::abs(m_data[ithread].v);
            if (tmp<val) val = tmp;
        }
        return val;
    }

};


#endif //M7_PRIVATESTORE_H
