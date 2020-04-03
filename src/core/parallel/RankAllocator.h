//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#include <src/core/table/Element.h>
#include "MPIWrapper.h"

template<typename T>
class RankAllocator {
    static_assert(std::is_base_of<Element, T>::value, "Rank allocation requires an Element-derived type");
    defs::inds m_data;
public:
    RankAllocator(size_t nblock) :m_data(nblock, 0ul) {
        size_t i = 0ul;
        for (auto &it:m_data) it = (i++)%mpi::nrank();
    }

    size_t get_rank(const T& key) const{
        return m_data[key.hash()%m_data.size()];
    }
};

#endif //M7_RANKALLOCATOR_H
