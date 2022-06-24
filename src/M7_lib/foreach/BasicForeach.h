//
// Created by Robert J. Anderson on 24/03/2022.
//

#ifndef M7_BASICFOREACH_H
#define M7_BASICFOREACH_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/Nd.h>
#include <M7_lib/util/Integer.h>
#include <M7_lib/util/Functor.h>
#include <M7_lib/util/Array.h>

class ExitLoop : public std::exception {
    virtual const char *what() const throw() {
        return "Loop body requested early termination of loop";
    }
};

namespace basic_foreach {
    /**
     * "compile time number of dimensions"
     */
    namespace ctnd {
        template<uint_t nind>
        using inds_t = defs::uinta_t<nind>;

        template<uint_t nind>
        struct Base {
        protected:
            static constexpr uint_t c_nind = nind;
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t<nind> m_value;
        public:
            /**
             * length of the index iterator
             */
            const uint_t m_niter;

            Base(uint_t niter) : m_niter(niter) {}
        };

        template<uint_t nind>
        struct Unrestricted : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_niter;
            const defs::uinta_t<nind> m_shape;
        private:
            /**
             * loop through the values of the ilevel element of m_value between 0 and the extent of that dimension
             * @tparam ilevel
             *  current level (position in the m_value array)
             */
            template<uint_t ilevel, typename fn_t>
            void level_loop(const fn_t &fn, tag::Int<ilevel>) {
                constexpr uint_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) level_loop(fn, tag::Int<ilevel + 1>());
            }

            /**
             * overload for when the last index has been reached
             */
            template<typename fn_t>
            void level_loop(const fn_t &fn, tag::Int<nind>) {
                constexpr uint_t iind = nind - 1;
                auto &ind = m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) fn(m_value);
            }

            /**
             * in the edge case that the nind is 0, do nothing
             */
            template<typename fn_t>
            void top_loop(const fn_t& /*fn*/, tag::Int<true>) {}

            /**
             * if nind is nonzero, start at the first index
             */
            template<typename fn_t>
            void top_loop(const fn_t& fn, tag::Int<false>) {
                level_loop(fn, tag::Int<1>());
            }

        public:

            static uint_t niter(const inds_t<nind> &shape){
                return nind ? nd::nelement(array::to_vector(shape)) : 0ul;
            }

            Unrestricted(const inds_t<nind> &shape) : Base<nind>(niter(shape)), m_shape(shape) {}

            Unrestricted(uint_t extent = 0ul) : Unrestricted(array::filled<uint_t, nind>(extent)) {}

            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor::assert_prototype<void(const inds_t<nind> &), fn_t>();
                top_loop(fn, tag::Int<nind == 0>());
            }
        };

        template<uint_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_niter;

            static uint_t niter(uint_t n) {
                if (!nind || !n) return 0ul;
                return integer::combinatorial(strict ? n : (n + nind) - 1, nind);
            }
        protected:
            uint_t m_n;

        private:

            template<uint_t ilevel, typename fn_t>
            void level_loop(const fn_t &fn, tag::Int<ilevel>) {
                constexpr uint_t iind = ascending ? (nind - ilevel) : (ilevel - 1);
                constexpr uint_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) level_loop(fn, tag::Int<ilevel + 1>());
            }

            template<typename fn_t>
            void level_loop(const fn_t &fn, tag::Int<nind>) {
                constexpr uint_t iind = ascending ? 0 : nind - 1;
                constexpr uint_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) fn(m_value);
            }

            template<typename fn_t>
            void top_loop(const fn_t& /*fn*/, tag::Int<true>) {}

            template<typename fn_t>
            void top_loop(const fn_t &fn, tag::Int<false>) {
                level_loop(fn, tag::Int<1>());
            }

        public:

            Ordered(uint_t n) : Base<nind>(niter(n)), m_n(n) {}

            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor::assert_prototype<void(const inds_t<nind> &), fn_t>();
                top_loop(fn, tag::Int<nind == 0>());
            }
        };
    }

    /**
     * "run time number of dimensions"
     */
    namespace rtnd {
        using inds_t = defs::uintv_t;

        struct Base {
            /**
             * length of the index iterator
             */
            const uint_t m_niter;
        protected:
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t m_value;
        public:
            /**
             * number of dimensions: length of the index array
             */
            const uint_t m_nind;

            Base(uint_t nind, uint_t niter) : m_niter(niter), m_value(nind, 0ul), m_nind(nind) {}
        };

        struct Unrestricted : Base {
            const inds_t m_shape;
        private:

            template<typename fn_t>
            void level_loop(const fn_t &fn, uint_t ilevel) {
                const auto &iind = ilevel - 1;
                auto &ind = m_value[iind];
                const auto &extent = m_shape[iind];
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(fn, ilevel + 1);
                else {
                    for (ind = 0ul; ind < extent; ++ind) fn(m_value);
                }
            }

        public:

            static uint_t niter(const inds_t& shape) {
                return shape.empty() ? 0ul : nd::nelement(shape);
            }

            static uint_t niter(uint_t nind, uint_t extent){
                return std::pow(extent, nind);
            }

            explicit Unrestricted(inds_t shape):
                Base(shape.size(), niter(shape)), m_shape(std::move(shape)) {}

            Unrestricted(uint_t nind, uint_t extent): Unrestricted(inds_t(nind, extent)) {}

            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor::assert_prototype<void(const inds_t &), fn_t>();
                if (!m_nind) return;
                level_loop(fn, 1);
            }
        };

        template<bool strict = true, bool ascending = true>
        struct Ordered : Base {

            static uint_t niter(uint_t n, uint_t r) {
                if (!r) return 0ul;
                return integer::combinatorial(strict ? n : (n + r) - 1, r);
            }

        private:
            const uint_t m_n;

            template<typename fn_t>
            void level_loop(const fn_t &fn, uint_t ilevel) {
                const uint_t iind = ascending ? (m_nind - ilevel) : (ilevel - 1);
                const uint_t iind_unrestrict = ascending ? m_nind - 1 : 0;
                auto &ind = m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(fn, ilevel + 1);
                else {
                    for (ind = 0ul; ind < extent; ++ind) fn(m_value);
                }
            }

        public:
            Ordered(uint_t n, uint_t r) : Base(r, niter(n, r)), m_n(n) {}

            template<typename fn_t>
            void loop(const fn_t& fn) {
                functor::assert_prototype<void(const inds_t &), fn_t>();
                if (!m_nind) return;
                level_loop(fn, 1);
            }
        };
    }
}


#endif //M7_BASICFOREACH_H
