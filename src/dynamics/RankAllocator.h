//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#include "src/hash/HashMap.h"

template<typename T>
class RankAllocator {
    typename HasherType<T> m_hasher;
    defs::inds m_data;
public:
    RankAllocator(size_t nblock) :m_data(n_block, 0ul) {
        size_t i = 0ul;
        for (auto &it:m_data) it = (i++)%mpi::nrank();
    }

    size_t get_rank(const T& key){
        return m_data[m_hasher(key)%m_data.size()];
    }
};


#endif //M7_RANKALLOCATOR_H
