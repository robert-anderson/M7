//
// Created by Robert J. Anderson on 19/11/2020.
//

#ifndef M7_REDUCTIONSYNDICATE_H
#define M7_REDUCTIONSYNDICATE_H

#if 0
#include "MPIWrapper.h"
#include "ReductionMemberBase.h"

struct ReductionSyndicateGroupBase {
    virtual void zero() = 0;

    virtual void sum(size_t iroot) = 0;

    virtual void all_sum() = 0;

    virtual void all_minloc() = 0;

    virtual void all_maxloc() = 0;
};

template<typename T>
class ReductionSyndicateGroup : public ReductionSyndicateGroupBase {

    std::vector<ReductionMemberBase<T> *> m_members;
    std::vector<T> m_local;
    std::vector<T> m_reduced;
    std::vector<size_t> m_rank_indices;
    std::vector<std::pair<T, size_t>> m_pairs;

    void split_pairs() {
        m_reduced.clear();
        m_rank_indices.clear();
        for (auto &pair: m_pairs) {
            m_reduced.push_back(pair.first);
            m_rank_indices.push_back(pair.second);
        }
    }

public:
    void add_member(ReductionMemberBase<T> *new_member) {
        m_members.push_back(new_member);
        m_local.resize(m_local.size() + new_member->m_nelement);
        m_reduced.resize(m_reduced.size() + new_member->m_nelement);
        m_rank_indices.resize(m_rank_indices.size() + new_member->m_nelement);
        m_pairs.resize(m_pairs.size() + new_member->m_nelement);
        ASSERT(m_pairs.size())
        size_t offset = 0ul;
        for (auto member: m_members) {
            member->update_data_ptrs(
                    (void *) (m_local.data() + offset),
                    (void *) (m_reduced.data() + offset),
                    (void *) (m_rank_indices.data() + offset));
            offset += member->m_nelement;
        }
    }

    void zero() override {
        std::memset(m_local.data(), 0, m_local.size() * sizeof(T));
        std::memset(m_reduced.data(), 0, m_reduced.size() * sizeof(T));
        std::memset(m_rank_indices.data(), 0, m_rank_indices.size() * sizeof(T));
    }

    void sum(size_t iroot) override {
        mpi::sum(m_local.data(), m_reduced.data(), m_local.size(), iroot);
    }

    void all_sum() override {
        mpi::all_sum(m_local.data(), m_reduced.data(), m_local.size());
    }

    void all_minloc() override {
        mpi::all_minloc(m_local, m_pairs);
        split_pairs();
    }

    void all_maxloc() override {
        mpi::all_maxloc(m_local, m_pairs);
        split_pairs();
    }

};

class ReductionSyndicate {

    std::array<std::unique_ptr<ReductionSyndicateGroupBase>, mpi_types.size()> m_groups;

public:
    template<typename T>
    void add_member(ReductionMemberBase<T> *new_member) {
        auto &group = m_groups[mpi_type_ind<T>()];
        if (!group.get()) {
            group = std::unique_ptr<ReductionSyndicateGroupBase>(new ReductionSyndicateGroup<T>);
        }
        static_cast<ReductionSyndicateGroup<T> *>(group.get())->add_member(new_member);
    }

    void zero() {
        for (auto &group: m_groups) if (group) group->zero();
    }

    void sum(size_t iroot = 0) {
        for (auto &group: m_groups) if (group) group->sum(iroot);
    }

    void all_sum() {
        for (auto &group: m_groups) if (group) group->all_sum();
    }

    void all_minloc() {
        for (auto &group: m_groups) if (group) group->all_minloc();
    }

    void all_maxloc() {
        for (auto &group: m_groups) if (group) group->all_maxloc();
    }

};


#endif //M7_REDUCTIONSYNDICATE_H
#endif //M7_REDUCTIONSYNDICATE_H