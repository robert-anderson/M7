//
// Created by Robert J. Anderson on 19/11/2020.
//

#ifndef M7_REDUCTIONMEMBERBASE_H
#define M7_REDUCTIONMEMBERBASE_H

#include <cstddef>

template<typename T>
struct ReductionMemberBase {
    size_t m_nelement;

    ReductionMemberBase(size_t nelement) : m_nelement(nelement) {}

    virtual void update_data_ptrs(void *local_ptr, void *reduced_ptr, void *rank_indices_ptr) = 0;
};


#endif //M7_REDUCTIONMEMBERBASE_H
