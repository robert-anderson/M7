//
// Created by Robert J. Anderson on 13/04/2021.
//

#ifndef M7_REDUCTION_H
#define M7_REDUCTION_H

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/util/SmartPtr.h"

template<typename T>
struct ReductionBase {
    const uint_t m_nelement;
    T *m_local_ptr = nullptr;
    T *m_reduced_ptr = nullptr;

    ReductionBase(uint_t nelement) : m_nelement(nelement) {}
};

template<typename T, uint_t nind>
struct NdReduction : ReductionBase<T> {
    typedef typename std::conditional<nind==0, T, buffered::Numbers<T, nind>>::type store_t;
    store_t m_local;
    store_t m_reduced;
    using ReductionBase<T>::m_local_ptr;
    using ReductionBase<T>::m_reduced_ptr;
    using ReductionBase<T>::m_nelement;

private:
    NdReduction(const NdFormat<nind>& /*format*/, tag::Int<true> /*zero_dims*/) :
            ReductionBase<T>(1) {
        m_local_ptr = &m_local;
        m_reduced_ptr = &m_reduced;
    }

    NdReduction(const NdFormat<nind>& format, tag::Int<false> /*zero_dims*/) :
            ReductionBase<T>(format.m_nelement), m_local(format.m_shape), m_reduced(format.m_shape) {
        m_local_ptr = reinterpret_cast<T *>(m_local.begin());
        m_reduced_ptr = reinterpret_cast<T *>(m_reduced.begin());
    }

public:
    NdReduction(const uinta_t<nind>& shape) :
            NdReduction({shape}, tag::Int<nind == 0>()){}

    NdReduction() : NdReduction({}, tag::Int<nind == 0>()){
        static_assert(!nind, "This ctor is only valid in the scalar case");
    }

    void all_sum() {
        mpi::all_sum(m_local_ptr, m_reduced_ptr, m_nelement);
    }

    void all_max() {
        mpi::all_max(m_local_ptr, m_reduced_ptr, m_nelement);
    }

    void all_min() {
        mpi::all_min(m_local_ptr, m_reduced_ptr, m_nelement);
    }
};

template<typename T>
using Reduction = NdReduction<T, 0>;


struct ReductionSyndicateGroupBase {

    virtual ~ReductionSyndicateGroupBase(){}

    virtual void zero_all_local() = 0;

    virtual void zero_all_reduced() = 0;

    virtual void zero_all() = 0;

    virtual void all_sum() = 0;

    virtual void all_max() = 0;

    virtual void all_min() = 0;
};

template<typename T>
struct ReductionSyndicateGroup : ReductionSyndicateGroupBase {
    std::vector<ReductionBase<T> *> m_members;
    std::vector<T> m_local_buffer;
    std::vector<T> m_reduced_buffer;

    ReductionSyndicateGroup() {}

    void add_member(ReductionBase<T>& member) {
        m_members.push_back(&member);
        m_local_buffer.resize(m_local_buffer.size() + member.m_nelement);
        m_reduced_buffer.resize(m_local_buffer.size());
    }

    void zero_all_local() override {
        for (ReductionBase<T> *member: m_members)
            std::memset(reinterpret_cast<char*>(member->m_local_ptr), 0, member->m_nelement * sizeof(T));
    }

    void zero_all_reduced() override {
        for (ReductionBase<T> *member: m_members)
            std::memset(reinterpret_cast<char*>(member->m_reduced_ptr), 0, member->m_nelement * sizeof(T));
    }

    void zero_all() override {
        zero_all_local();
        zero_all_reduced();
    }

private:
    void collect() {
        auto dst = m_local_buffer.data();
        for (ReductionBase<T> *member: m_members) {
            const uint_t nelement = member->m_nelement;
            std::memcpy(dst, member->m_local_ptr, nelement * sizeof(T));
            dst += nelement;
        }
    }

    void disperse() {
        auto src = m_reduced_buffer.data();
        for (ReductionBase<T> *member: m_members) {
            const uint_t nelement = member->m_nelement;
            std::memcpy(member->m_reduced_ptr, src, nelement * sizeof(T));
            src += nelement;
        }
    }

public:
    void all_sum() override {
        collect();
        mpi::all_sum(m_local_buffer.data(), m_reduced_buffer.data(), m_local_buffer.size());
        disperse();
    }

    void all_max() override {
        collect();
        mpi::all_max(m_local_buffer.data(), m_reduced_buffer.data(), m_local_buffer.size());
        disperse();
    }

    void all_min() override {
        collect();
        mpi::all_min(m_local_buffer.data(), m_reduced_buffer.data(), m_local_buffer.size());
        disperse();
    }
};

struct ReductionSyndicate {
    std::array<std::unique_ptr<ReductionSyndicateGroupBase>, mpi_types.size()> m_groups;

    template<typename T>
    void add_member(ReductionBase<T>& member) {
        auto igroup = mpi_type_ind<T>();
        if (!m_groups[igroup])
            m_groups[igroup] = smart_ptr::make_unique<ReductionSyndicateGroup<T>>();
        dynamic_cast<ReductionSyndicateGroup<T> *>(m_groups[igroup].get())->add_member(member);
    }

    void add_members() {}

    template<typename T, typename ...Args>
    void add_members(ReductionBase<T>& first, Args& ... rest) {
        add_member(first);
        add_members(std::forward<Args&>(rest)...);
    }

    void zero_all_local() {
        for (auto& group: m_groups) if (group) group->zero_all_local();
    }

    void zero_all_reduced() {
        for (auto& group: m_groups) if (group) group->zero_all_reduced();
    }

    void zero_all() {
        for (auto& group: m_groups) if (group) group->zero_all();
    }

    void all_sum() {
        for (auto& group: m_groups) if (group) group->all_sum();
    }

    void all_max() {
        for (auto& group: m_groups) if (group) group->all_max();
    }

    void all_min() {
        for (auto& group: m_groups) if (group) group->all_min();
    }
};

#endif //M7_REDUCTION_H
