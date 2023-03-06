//
// Created by Robert J. Anderson on 13/04/2021.
//

#ifndef M7_REDUCTION_H
#define M7_REDUCTION_H

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/util/Pointer.h"

namespace reduction {

    template<typename T>
    struct ReductionBase {
        const uint_t m_nelement;
        T *m_local_ptr = nullptr;
        T *m_reduced_ptr = nullptr;

        ReductionBase(uint_t nelement) : m_nelement(nelement) {}
    };

    template<typename T, uint_t nind>
    struct LocalReducedPair : ReductionBase<T> {
        typedef typename std::conditional<nind==0, T, buffered::Numbers<T, nind>>::type store_t;
        store_t m_local;
        store_t m_reduced;
        using ReductionBase<T>::m_local_ptr;
        using ReductionBase<T>::m_reduced_ptr;
        using ReductionBase<T>::m_nelement;

    private:
        LocalReducedPair(const NdFormat<nind>& /*format*/, tag::Int<true> /*zero_dims*/) :
                ReductionBase<T>(1) {
            m_local_ptr = &m_local;
            m_reduced_ptr = &m_reduced;
            m_local = 0;
            m_reduced = 0;
        }

        LocalReducedPair(const NdFormat<nind>& format, tag::Int<false> /*zero_dims*/) :
                ReductionBase<T>(format.m_nelement), m_local(format.m_shape), m_reduced(format.m_shape) {
            m_local_ptr = reinterpret_cast<T *>(m_local.begin());
            m_reduced_ptr = reinterpret_cast<T *>(m_reduced.begin());
            m_local = 0;
            m_reduced = 0;
        }

    public:
        LocalReducedPair(const uinta_t<nind>& shape) :
                LocalReducedPair({shape}, tag::Int<nind == 0>()){}

        LocalReducedPair() : LocalReducedPair({}, tag::Int<nind == 0>()){
            static_assert(!nind, "This ctor is only valid in the scalar case");
        }
    };


    template<typename T, uint_t nind>
    struct NdArray : LocalReducedPair<T, nind> {
        using ReductionBase<T>::m_local_ptr;
        using ReductionBase<T>::m_reduced_ptr;
        using ReductionBase<T>::m_nelement;

        NdArray(const uinta_t<nind>& shape) : LocalReducedPair<T, nind>({shape}){}
        NdArray() : LocalReducedPair<T, 0ul>(){}

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
    using Scalar = NdArray<T, 0>;


    namespace cyclic {

        struct CyclicBase {
            virtual void after_delta_reduce() {}
        };

        template<typename T>
        struct Cyclic : CyclicBase {
            typedef dtype::make_signed_if_integral<T> delta_t;
            ReductionBase<delta_t>* m_delta_base;
            ReductionBase<T>* m_total_base;
            ReductionBase<delta_t>* m_prev_delta_base;
            ReductionBase<T>* m_prev_total_base;

            Cyclic(ReductionBase<delta_t>* delta_base, ReductionBase<T>* total_base,
                   ReductionBase<delta_t>* prev_delta_base, ReductionBase<T>* prev_total_base):
                m_delta_base(delta_base), m_total_base(total_base),
                m_prev_delta_base(prev_delta_base), m_prev_total_base(prev_total_base){}

            void after_delta_reduce() override {
                const auto nelement = m_total_base->m_nelement;
                // update previous totals
                {
                    auto src = m_total_base->m_local_ptr;
                    auto dst = m_prev_total_base->m_local_ptr;
                    std::copy(src, src + nelement, dst);
                }
                {
                    auto src = m_total_base->m_reduced_ptr;
                    auto dst = m_prev_total_base->m_reduced_ptr;
                    std::copy(src, src + nelement, dst);
                }
                // update previous deltas
                {
                    auto src = m_delta_base->m_local_ptr;
                    auto dst = m_prev_delta_base->m_local_ptr;
                    std::copy(src, src + nelement, dst);
                }
                {
                    auto src = m_delta_base->m_reduced_ptr;
                    auto dst = m_prev_delta_base->m_reduced_ptr;
                    std::copy(src, src + nelement, dst);
                }
                // update totals
                for (uint_t ielement=0ul; ielement < nelement; ++ielement) {
                    m_total_base->m_local_ptr[ielement] += m_delta_base->m_local_ptr[ielement];
                    m_total_base->m_reduced_ptr[ielement] += m_delta_base->m_reduced_ptr[ielement];
                }
                // re-zero deltas
                {
                    auto ptr = m_delta_base->m_local_ptr;
                    std::memset(ptr, 0, nelement * sizeof(delta_t));
                }
                {
                    auto ptr = m_delta_base->m_reduced_ptr;
                    std::memset(ptr, 0, nelement * sizeof(delta_t));
                }
            }
        };

        template<typename T, uint_t nind>
        struct NdArray : Cyclic<T> {
            using delta_t = typename Cyclic<T>::delta_t;
        private:
            reduction::NdArray<delta_t, nind> m_delta;
            // total should never be directly modified - the point of this class is to modify by deltas and accumulate
            reduction::LocalReducedPair<T, nind> m_total;
            reduction::LocalReducedPair<delta_t, nind> m_prev_delta;
            reduction::LocalReducedPair<T, nind> m_prev_total;
            // don't expose these pointers to the base classes - only needed for adding as member of a syndicate
            using Cyclic<T>::m_delta_base;
            using Cyclic<T>::m_total_base;
            using Cyclic<T>::after_delta_reduce;

        public:
            /**
             * @return
             *  modifiable ref to the local delta array
             */
            typename reduction::LocalReducedPair<delta_t, nind>::store_t& delta() {
                return m_delta.m_local;
            }

            /**
             * @return
             *  const ref to the reduced total
             */
            const typename reduction::LocalReducedPair<T, nind>::store_t& total() const {
                return m_total.m_reduced;
            }

            const typename reduction::LocalReducedPair<delta_t, nind>::store_t& prev_delta() const {
                return m_prev_delta.m_reduced;
            }

            const typename reduction::LocalReducedPair<T, nind>::store_t& prev_total() const {
                return m_prev_total.m_reduced;
            }


            NdArray(const uinta_t<nind>& shape) :
                Cyclic<T>(&m_delta, &m_total, &m_prev_delta, &m_prev_total),
                m_delta(shape), m_total(shape), m_prev_delta(shape), m_prev_total(shape){}

            NdArray(): NdArray(uinta_t<0ul>{}){}

            // cyclic only needs to support accumulations i.e. summation
            void all_sum() {
                m_delta.all_sum();
                after_delta_reduce();
            }
        };

        template<typename T>
        using Scalar = NdArray<T, 0>;
    }

    struct SyndicateGroupBase {

        virtual ~SyndicateGroupBase(){}

        virtual void zero_all_local() = 0;

        virtual void zero_all_reduced() = 0;

        virtual void zero_all() = 0;

        virtual void all_sum() = 0;

        virtual void all_max() = 0;

        virtual void all_min() = 0;
    };

    template<typename T>
    struct SyndicateGroup : reduction::SyndicateGroupBase {
        v_t<ReductionBase<T> *> m_members;
        v_t<T> m_local_buffer;
        v_t<T> m_reduced_buffer;

        SyndicateGroup() {}

        void add_member(ReductionBase<T>& member) {
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
                memset(reinterpret_cast<char*>(member->m_reduced_ptr), 0, member->m_nelement * sizeof(T));
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
                memcpy(dst, member->m_local_ptr, nelement * sizeof(T));
                dst += nelement;
            }
        }

        void disperse() {
            auto src = m_reduced_buffer.data();
            for (ReductionBase<T> *member: m_members) {
                const uint_t nelement = member->m_nelement;
                memcpy(member->m_reduced_ptr, src, nelement * sizeof(T));
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

    struct Syndicate {
        std::array<std::unique_ptr<SyndicateGroupBase>, mpi_types.size()> m_groups;
        v_t<cyclic::CyclicBase*> m_cyclic_members;

    private:
        /**
         * @return
         *  syndicate group of the required type. if it does not already exist, it will be created first
         */
        template<typename T>
        SyndicateGroup<T>* get_group() {
            auto igroup = mpi_type_ind<T>();
            if (!m_groups[igroup]) m_groups[igroup] = ptr::smart::make_unique<SyndicateGroup<T>>();
            return dynamic_cast<SyndicateGroup<T> *>(m_groups[igroup].get());
        }

    public:

        template<typename T>
        void add_member(ReductionBase<T>& member) {
            get_group<T>()->add_member(member);
        }

        template<typename T>
        void add_member(cyclic::Cyclic<T>& member) {
            add_member(*member.m_delta_base);
            m_cyclic_members.push_back(&member);
        }

        void add_members() {}

        template<typename T, typename ...Args>
        void add_members(ReductionBase<T>& first, Args& ... rest) {
            add_member(first);
            add_members(std::forward<Args&>(rest)...);
        }

        template<typename T, typename ...Args>
        void add_members(cyclic::Cyclic<T>& first, Args& ... rest) {
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
            for (auto cyclic_member: m_cyclic_members) cyclic_member->after_delta_reduce();
        }

        void all_max() {
            for (auto& group: m_groups) if (group) group->all_max();
            for (auto cyclic_member: m_cyclic_members) cyclic_member->after_delta_reduce();
        }

        void all_min() {
            for (auto& group: m_groups) if (group) group->all_min();
            for (auto cyclic_member: m_cyclic_members) cyclic_member->after_delta_reduce();
        }
    };
}


#endif //M7_REDUCTION_H
