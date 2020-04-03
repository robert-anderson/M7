//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#if 0
#include "src/hash/HashMap.h"

template<typename T>
class RankAllocator {
    typename HasherType<T>::type m_hasher;
    defs::inds m_data;
public:
    RankAllocator(size_t nblock) :m_data(nblock, 0ul) {
        size_t i = 0ul;
        for (auto &it:m_data) it = (i++)%mpi::nrank();
    }

    size_t get_rank(const T& key) const{
        //return m_data[m_hasher(key)%m_data.size()];
    }
};


#endif //M7_RANKALLOCATOR_H
#endif //M7_RANKALLOCATOR_H
