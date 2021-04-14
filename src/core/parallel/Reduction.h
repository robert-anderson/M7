//
// Created by rja on 13/04/2021.
//

#ifndef M7_REDUCTION_H
#define M7_REDUCTION_H

template<typename T>
struct ReductionBase {
    const size_t m_nelement;
    T *m_local_ptr = nullptr;
    T *m_reduced_ptr = nullptr;

    ReductionBase(size_t nelement) : m_nelement(nelement) {}
};

template<typename T, size_t nind>
struct Reduction : ReductionBase<T> {
    buffered::Numbers<T, nind> m_local;
    buffered::Numbers<T, nind> m_reduced;
    using ReductionBase<T>::m_local_ptr;
    using ReductionBase<T>::m_reduced_ptr;
    using ReductionBase<T>::m_nelement;


    Reduction(const NdFormat<nind> &format) :
            ReductionBase<T>(format.nelement()), m_local(format.shape()), m_reduced(format.shape()) {
        m_local_ptr = reinterpret_cast<T *>(m_local.begin());
        m_reduced_ptr = reinterpret_cast<T *>(m_reduced.begin());
    }

    void all_sum() {
        mpi::all_sum(m_local_ptr, m_reduced_ptr, m_nelement);
    }
};

struct ReductionSyndicateGroupBase {
    virtual void zero_all_local() = 0;

    virtual void zero_all_reduced() = 0;

    virtual void zero_all() = 0;

    virtual void all_sum() = 0;
};

template<typename T>
struct ReductionSyndicateGroup : ReductionSyndicateGroupBase {
    std::vector<ReductionBase<T> *> m_members;
    std::vector<T> m_local_buffer;
    std::vector<T> m_reduced_buffer;

    ReductionSyndicateGroup() {}

    void add_member(ReductionBase<T> &member) {
        m_members.push_back(&member);
        m_local_buffer.resize(m_local_buffer.size() + member.m_nelement);
        m_reduced_buffer.resize(m_local_buffer.size());
    }

    void zero_all_local() override {
        for (ReductionBase<T> *member: m_members)
            std::memset(member->m_local_ptr, 0, member->m_nelement * sizeof(T));
    }

    void zero_all_reduced() override {
        for (ReductionBase<T> *member: m_members)
            std::memset(member->m_reduced_ptr, 0, member->m_nelement * sizeof(T));
    }

    void zero_all() override {
        zero_all_local();
        zero_all_reduced();
    }

private:
    void collect() {
        auto dst = m_local_buffer.data();
        for (ReductionBase<T> *member: m_members) {
            const size_t nelement = member->m_nelement;
            std::memcpy(dst, member->m_local_ptr, nelement * sizeof(T));
            dst += nelement;
        }
    }

    void disperse() {
        auto src = m_reduced_buffer.data();
        for (ReductionBase<T> *member: m_members) {
            const size_t nelement = member->m_nelement;
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
};

struct ReductionSyndicate {
    std::array<std::unique_ptr<ReductionSyndicateGroupBase>, mpi_types.size()> m_groups;

    template<typename T>
    void add_member(ReductionBase<T> &member) {
        auto igroup = mpi_type_ind<T>();
        if (!m_groups[igroup])
            m_groups[igroup] = std::unique_ptr<ReductionSyndicateGroup<T>>(new ReductionSyndicateGroup<T>());
        dynamic_cast<ReductionSyndicateGroup<T> *>(m_groups[igroup].get())->add_member(member);
    }

    void add_members() {}

    template<typename T, typename ...Args>
    void add_members(ReductionBase<T> &first, Args &... rest) {
        add_member(first);
        add_members(std::forward<Args &>(rest)...);
    }

    void zero_all_local() {
        for (auto &group: m_groups) if (group) group->zero_all_local();
    }

    void zero_all_reduced() {
        for (auto &group: m_groups) if (group) group->zero_all_reduced();
    }

    void zero_all() {
        for (auto &group: m_groups) if (group) group->zero_all();
    }

    void all_sum() {
        for (auto &group: m_groups) if (group) group->all_sum();
    }
};

#endif //M7_REDUCTION_H