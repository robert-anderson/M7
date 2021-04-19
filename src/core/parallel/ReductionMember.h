//
// Created by RJA on 19/11/2020.
//

#ifndef M7_REDUCTIONMEMBER_H
#define M7_REDUCTIONMEMBER_H

#if 0
#include "src/core/nd/NdAccessor.h"
#include "ReductionSyndicate.h"


template<typename T, size_t nind=0>
struct ReductionMember : FormattedNdAccessor<T, nind>, ReductionMemberBase<T> {
    FormattedNdAccessor<T, nind> m_reduced;
    FormattedNdAccessor<size_t, nind> m_rank_indices;

    ReductionMember(ReductionSyndicate& syndicate, NdFormat<nind> format) :
            FormattedNdAccessor<T, nind>(nullptr, format),
            ReductionMemberBase<T>(format.nelement()),
            m_reduced(nullptr, format), m_rank_indices(nullptr, format) {
        syndicate.add_member(this);
    }

    void update_data_ptrs(void *local_ptr, void *reduced_ptr, void* rank_indices_ptr) override {
        NdAccessor<T, nind>::m_data = (T*)local_ptr;
        m_reduced.m_data = (T*)reduced_ptr;
        m_rank_indices.m_data = (size_t *)rank_indices_ptr;
    }

    template<typename ...Args>
    const T &reduced(Args... inds) const {
        return m_reduced(inds...);
    }

    template<typename ...Args>
    const size_t &rank_index(Args... inds) const {
        return m_rank_indices(inds...);
    }
};

template<typename T, size_t nind=0>
struct SingleReducible : ReductionSyndicate, ReductionMember<T, nind> {
    SingleReducible(NdFormat<nind> format) : ReductionMember<T, nind>(*this, format){}
};


#endif //M7_REDUCTIONMEMBER_H
#endif //M7_REDUCTIONMEMBER_H
