//
// Created by rja on 24/03/2022.
//

#ifndef M7_BASICFOREACH_H
#define M7_BASICFOREACH_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/utils.h>
#include <M7_lib/nd/NdFormat.h>
#include <M7_lib/nd/NdFormatD.h>

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
        template<size_t nind>
        using inds_t = std::array<size_t, nind>;

        template<size_t nind>
        struct Base {
        protected:
            static constexpr size_t c_nind = nind;
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t<nind> m_value;
        public:
            /**
             * length of the index iterator
             */
            const size_t m_niter;


            Base(size_t niter) : m_niter(niter) {}
        };

        template<size_t nind>
        struct Unrestricted : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_niter;
            const inds_t<nind> m_shape;
        private:
            /**
             * loop through the values of the ilevel element of m_value between 0 and the extent of that dimension
             * @tparam ilevel
             *  current level (position in the m_value array)
             */
            template<size_t ilevel, typename fn_t>
            void level_loop(const fn_t &fn, tags::Ind<ilevel>) {
                constexpr size_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) level_loop(fn, tags::Ind<ilevel + 1>());
            }

            /**
             * overload for when the last index has been reached
             */
            template<typename fn_t>
            void level_loop(const fn_t &fn, tags::Ind<nind>) {
                constexpr size_t iind = nind - 1;
                auto &ind = m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) fn(m_value);
            }

            /**
             * in the edge case that the nind is 0, do nothing
             */
            template<typename fn_t>
            void top_loop(const fn_t &fn, tags::Bool<true>) {}

            /**
             * if nind is nonzero, start at the first index
             */
            template<typename fn_t>
            void top_loop(const fn_t &fn, tags::Bool<false>) {
                level_loop(fn, tags::Ind<1>());
            }

        public:

            Unrestricted(const inds_t<nind> &shape) : Base<nind>(NdFormat<nind>(shape).m_nelement), m_shape(shape) {}

            Unrestricted(size_t extent = 0ul) : Unrestricted(array_utils::filled<size_t, nind>(extent)) {}

            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const inds_t<nind> &), fn_t>();
                top_loop(fn, tags::Bool<nind == 0>());
            }
        };

        template<size_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_niter;
        protected:
            size_t m_n;

        private:
            static size_t make_nterm(size_t n) {
                if (!nind) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + nind) - 1, nind);
            }

            template<size_t ilevel, typename fn_t>
            void level_loop(const fn_t &fn, tags::Ind<ilevel>) {
                constexpr size_t iind = ascending ? (nind - ilevel) : (ilevel - 1);
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) level_loop(fn, tags::Ind<ilevel + 1>());
            }

            template<typename fn_t>
            void level_loop(const fn_t &fn, tags::Ind<nind>) {
                constexpr size_t iind = ascending ? 0 : nind - 1;
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) fn(m_value);
            }

            template<typename fn_t>
            void top_loop(const fn_t &fn, tags::Bool<true>) {}

            template<typename fn_t>
            void top_loop(const fn_t &fn, tags::Bool<false>) {
                level_loop(fn, tags::Ind<1>());
            }

        public:

            Ordered(size_t n) : Base<nind>(make_nterm(n)), m_n(n) {}

            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const inds_t<nind> &), fn_t>();
                top_loop(fn, tags::Bool<nind == 0>());
            }
        };
    }

    /**
     * "run time number of dimensions"
     */
    namespace rtnd {
        using inds_t = std::vector<size_t>;

        struct Base {
            /**
             * length of the index iterator
             */
            const size_t m_niter;
        protected:
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t m_value;
        public:
            /**
             * number of dimensions: length of the index array
             */
            const size_t m_nind;

            Base(size_t nind, size_t niter) : m_niter(niter), m_value(nind, 0ul), m_nind(nind) {}
        };

        struct Unrestricted : Base {
            const inds_t m_shape;
        private:

            template<typename fn_t>
            void level_loop(const fn_t &fn, size_t ilevel) {
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
            explicit Unrestricted(inds_t shape):
                Base(shape.size(), NdFormatD(shape).m_nelement), m_shape(std::move(shape)) {}

            Unrestricted(size_t nind, size_t extent):
                    Unrestricted(std::vector<size_t>(nind, extent)) {}

            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const inds_t &), fn_t>();
                if (!m_nind) return;
                level_loop(fn, 1);
            }
        };

        template<bool strict = true, bool ascending = true>
        struct Ordered : Base {
        private:
            const size_t m_n;

            static size_t make_nterm(size_t n, size_t r) {
                if (!r) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + r) - 1, r);
            }

            template<typename fn_t>
            void level_loop(const fn_t &fn, size_t ilevel) {
                const size_t iind = ascending ? (m_nind - ilevel) : (ilevel - 1);
                const size_t iind_unrestrict = ascending ? m_nind - 1 : 0;
                auto &ind = m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(fn, ilevel + 1);
                else {
                    for (ind = 0ul; ind < extent; ++ind) fn(m_value);
                }
            }

        public:
            Ordered(size_t n, size_t r) : Base(r, make_nterm(n, r)), m_n(n) {}

            template<typename fn_t>
            void loop(const fn_t& fn) {
                functor_utils::assert_prototype<void(const inds_t &), fn_t>();
                if (!m_nind) return;
                level_loop(fn, 1);
            }
        };
    }
}


#endif //M7_BASICFOREACH_H