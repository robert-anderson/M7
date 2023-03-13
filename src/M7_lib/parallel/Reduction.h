//
// Created by Robert J. Anderson on 13/04/2021.
//

#ifndef M7_REDUCTION_H
#define M7_REDUCTION_H

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/util/Pointer.h"

namespace reduction {

    struct Base {
        FieldBase& m_local_base;
        FieldBase& m_reduced_base;
        const uint_t m_itype;
        Base(FieldBase& local_base, FieldBase& reduced_base, uint_t itype):
            m_local_base(local_base), m_reduced_base(reduced_base), m_itype(itype){}

        virtual ~Base() = default;

        virtual void post_reduce(){}

    private:
        uint_t m_global_reduced_offset = ~0ul;

    public:
        void to_global_send() {
            auto& send_buf = g_send_reduction_buffers[m_itype];
            m_global_reduced_offset = send_buf.size();
            auto ptr = m_local_base.cbegin();
            send_buf.insert(send_buf.cend(), ptr, ptr + m_local_base.m_size);
            auto& recv_buf = g_recv_reduction_buffers[m_itype];
            if (recv_buf.size() < send_buf.size()) recv_buf.resize(send_buf.size());
        }

        void from_global_recv() {
            DEBUG_ASSERT_NE(m_global_reduced_offset, ~0ul, "global offset should have been set in a prior to_global_send call");
            auto& buf = g_recv_reduction_buffers[m_itype];
            const auto ptr = m_reduced_base.begin();
            std::memcpy(ptr, buf.data()+m_global_reduced_offset, m_reduced_base.m_size);
            m_global_reduced_offset = ~0ul;
        }
    };

    template<typename T, uint_t nind>
    struct NdArray : public Base {
        /*
         * buffered::Number class has additional methods for conversion
         */
        typedef typename std::conditional<nind==0, buffered::Number<T>, buffered::Numbers<T, nind>>::type store_t;
        NdFormat<nind> m_format;
        store_t m_local;
        store_t m_reduced;
        NdArray(const uinta_t<nind>& shape) :
            Base(m_local, m_reduced, mpi_type_ind<T>()), m_format(shape), m_local(m_format), m_reduced(m_format){}

        NdArray& operator=(const NdArray<T, nind>& other) {
            m_local = other.m_local;
            m_reduced = other.m_reduced;
            return *this;
        }

        NdArray(const NdArray<T, nind>& other) : NdArray(other.m_format.m_shape){
            *this = other;
        }

        // only available for scalar case
        NdArray() : Base(m_local, m_reduced, mpi_type_ind<T>()){}
    private:
        void all_reduce(MpiOp op) {
            mpi::all_reduce(m_local.ctbegin(), m_reduced.tbegin(), op, m_local.m_nelement);
            post_reduce();
        }

    public:
        void all_sum() {all_reduce(MpiSum);}
        void all_max() {all_reduce(MpiMax);}
        void all_min() {all_reduce(MpiMin);}
    };

    template<typename T>
    struct Scalar : NdArray<T, 0> {
        Scalar() : NdArray<T, 0>() {}

        Scalar& operator=(const Scalar& other) {
            NdArray<T, 0>::operator=(other);
            return *this;
        }

        Scalar(const Scalar<T>& other) : NdArray<T, 0>() {
            *this = other;
        }
    };


    namespace cyclic {

        template<typename delta_t, uint_t nind, bool total_signed = std::is_signed<delta_t>::value>
        struct NdArray : reduction::NdArray<delta_t, nind> {
            typedef dtype::set_signedness_if_integral_t<delta_t, total_signed> total_t;

        protected:
            using reduction::NdArray<delta_t, nind>::m_local;
            using reduction::NdArray<delta_t, nind>::m_reduced;

            reduction::NdArray<total_t, nind> m_total;
            reduction::NdArray<delta_t, nind> m_prev_delta;
            reduction::NdArray<total_t, nind> m_prev_total;


            // only available for scalar case
            NdArray(): reduction::NdArray<delta_t, 0ul>() {}
        public:
            NdArray(const uinta_t<nind>& shape):
                reduction::NdArray<delta_t, nind>(shape), m_total(shape), m_prev_delta(shape), m_prev_total(shape){}


            void post_reduce() override {
                const auto nelement = m_local.m_nelement;
                // update previous deltas
                m_prev_delta = *this;
                // update previous totals
                m_prev_total = m_total;
                // update totals
                for (uint_t ielement=0ul; ielement < nelement; ++ielement) {
                    m_total.m_local.tbegin()[ielement] += m_local.ctbegin()[ielement];
                    m_total.m_reduced.tbegin()[ielement] += m_reduced.ctbegin()[ielement];
                }
                // re-zero deltas
                m_local = 0;
                m_reduced = 0;
            }

            const typename reduction::NdArray<total_t, nind>::store_t& total() const {
                return m_total.m_reduced;
            }

            const typename reduction::NdArray<total_t, nind>::store_t& prev_total() const {
                return m_prev_total.m_reduced;
            }

            typename reduction::NdArray<delta_t, nind>::store_t& delta() {
                return m_local;
            }

            const reduction::NdArray<delta_t, nind>& prev_delta() const {
                return m_prev_delta;
            }
        };

        template<typename delta_t, bool total_signed = std::is_signed<delta_t>::value>
        struct Scalar : NdArray<delta_t, 0, total_signed> {
            Scalar() : NdArray<delta_t, 0, total_signed>() {}

            Scalar& operator=(const Scalar& other) {
                NdArray<delta_t, 0, total_signed>::operator=(other);
                return *this;
            }

            Scalar(const Scalar<delta_t, total_signed>& other) : NdArray<delta_t, 0, total_signed>() {
                *this = other;
            }
        };
    }

    namespace {
        template<typename T>
        void all_reduce_one_type(MpiOp op) {
            const auto& send = g_send_reduction_buffers[mpi_type_ind<T>()];
            auto& recv = g_recv_reduction_buffers[mpi_type_ind<T>()];
            DEBUG_ASSERT_EQ(send.size(), recv.size(), "send and recv buffers for the same type should be identical");
            if (!send.size()) return;
            auto send_ptr = reinterpret_cast<const T*>(send.data());
            auto recv_ptr = reinterpret_cast<T*>(recv.data());
            mpi::all_reduce(send_ptr, recv_ptr, op, send.size() / sizeof(T));
        }

        void all_reduce(MpiOp op) {
            all_reduce_one_type<char>(op);
            all_reduce_one_type<short int>(op);
            all_reduce_one_type<int>(op);
            all_reduce_one_type<long int>(op);
            all_reduce_one_type<long long int>(op);
            all_reduce_one_type<unsigned char>(op);
            all_reduce_one_type<unsigned short int>(op);
            all_reduce_one_type<unsigned int>(op);
            all_reduce_one_type<unsigned long int>(op);
            all_reduce_one_type<unsigned long long int>(op);
            all_reduce_one_type<float>(op);
            all_reduce_one_type<double>(op);
            all_reduce_one_type<long double>(op);
            all_reduce_one_type<std::complex<float>>(op);
            all_reduce_one_type<std::complex<double>>(op);
            all_reduce_one_type<std::complex<long double>>(op);
            all_reduce_one_type<bool>(op);
        }

        void all_reduce(v_t<Base*>& members, MpiOp op) {
            for (auto& member: members) member->to_global_send();
            all_reduce(op);
            for (auto& member: members) member->from_global_recv();
            for (auto& send: g_send_reduction_buffers) send.clear();
            for (auto& recv: g_recv_reduction_buffers) recv.assign(recv.size(), 0);
        }
    }

    static void all_sum(v_t<Base*>& members) {all_reduce(members, MpiSum);}
    static void all_max(v_t<Base*>& members) {all_reduce(members, MpiMax);}
    static void all_min(v_t<Base*>& members) {all_reduce(members, MpiMin);}

    static void zero_local(v_t<Base*>& members) {
        for (auto& member: members) member->m_local_base.zero();
    }

    static void all_reduce(v_t<Base*>& members, MpiOp op) {
        for (auto& member: members) member->to_global_send();
        all_reduce(op);
        for (auto& member: members) member->from_global_recv();
        for (auto& send: g_send_reduction_buffers) send.clear();
        for (auto& recv: g_recv_reduction_buffers) recv.assign(recv.size(), 0);
    }
}

#endif //M7_REDUCTION_H