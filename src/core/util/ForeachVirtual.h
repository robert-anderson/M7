//
// Created by anderson on 2/9/22.
//

#ifndef M7_FOREACHVIRTUAL_H
#define M7_FOREACHVIRTUAL_H

#include <src/core/parallel/MPIAssert.h>
#include "utils.h"

namespace foreach_virtual {

    class ExitLoop : public std::exception {
        virtual const char *what() const throw() {
            return "Loop body requested early termination of loop";
        }
    };

    /**
     * "compile time number of dimensions"
     */
    namespace ctnd {

        template<size_t nind>
        using inds_t = std::array<size_t, nind>;

        template<size_t nind>
        struct Base {
        protected:
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t<nind> m_value;
            /**
             * integer index of the current iteration (body call)
             */
            size_t m_iiter;
        public:
            /**
             * getter for the current index value so that body implementations make immutable reference to the indices
             */
            const inds_t<nind> &value() const { return m_value; }

            /**
             * getter for the current iteration index
             */
            const size_t &iiter() const { return m_iiter; }

            /**
             * length of the index iterator
             */
            const size_t m_niter;

            /**
             * @return
             *  sum over current indices
             */
            size_t sum() const {
                return std::accumulate(m_value.cbegin(), m_value.cend(), 0);
            }

            void set_null_value() {
                m_value.fill(~0ul);
            }

            Base(size_t nterm) : m_niter(nterm) {}

            /**
             * function to be called each time a new set of indices is formed
             */
            virtual void body() = 0;

        protected:
            /**
             * function defining the algorithm by which sets of indices are generated
             */
            virtual void throwing_loop() = 0;

        public:
            void loop() {
                if (nind) {
                    m_iiter = ~0ul;
                    try {
                        throwing_loop();
                        // if loop completes without early termination, increment to final value
                        ++m_iiter;
                        DEBUG_ASSERT_EQ(m_iiter, m_niter, "loop completed after incorrect number of iterations");
                    }
                    catch (const ExitLoop &) {}
                } else m_iiter = 0ul;
                // leave value "undefined"
                set_null_value();
            }
        };

        template<size_t nind>
        struct Unrestricted : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_iiter;
            using Base<nind>::body;
            const inds_t<nind> m_shape;
        private:
            /**
             * @param shape
             *  extents of each index
             * @return
             *  expected number of calls to body() without early termination
             */
            static size_t nterm(const inds_t<nind> &shape) {
                if (!nind) return 0;
                size_t n = 1;
                for (size_t i = 0ul; i != nind; ++i) n *= shape[i];
                return n;
            }

            /**
             * loop through the values of the ilevel element of m_value between 0 and the extent of that dimension
             * @tparam ilevel
             *  current level (position in the m_value array)
             */
            template<size_t ilevel>
            void level_loop(tags::Ind<ilevel>) {
                constexpr size_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    try { level_loop(tags::Ind<ilevel + 1>()); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            /**
             * overload for when the last index has been reached
             */
            void level_loop(tags::Ind<nind>) {
                constexpr size_t iind = nind - 1;
                auto &ind = m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    ++m_iiter;
                    try { body(); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            /**
             * in the edge case that the nind is 0, do nothing but set the iteration counter
             */
            void top_loop(tags::Bool<true>) {
                m_iiter = 0ul;
            }

            /**
             * if nind is nonzero, start at the first index
             */
            void top_loop(tags::Bool<false>) {
                try { level_loop(tags::Ind<1>()); }
                catch (const ExitLoop &ex) {throw ex;}
            }

        public:
            Unrestricted(const inds_t<nind> &shape) : Base<nind>(nterm(shape)), m_shape(shape) {}

            Unrestricted(const size_t &extent) : Unrestricted(array_utils::filled<size_t, nind>(extent)) {}

        protected:
            void throwing_loop() override {
                top_loop(tags::Bool<nind == 0>());
            }
        };

        template<size_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_iiter;
            using Base<nind>::body;
            using inds_t = std::array<size_t, nind>;
            const size_t m_n;

        private:
            static size_t nterm(size_t n) {
                if (!nind) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + nind) - 1, nind);
            }

            template<size_t ilevel>
            void level_loop(tags::Ind<ilevel>) {
                constexpr size_t iind = ascending ? (nind - ilevel) : (ilevel - 1);
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    try { level_loop(tags::Ind<ilevel + 1>()); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            void level_loop(tags::Ind<nind>) {
                constexpr size_t iind = ascending ? 0 : nind - 1;
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    ++m_iiter;
                    try { body(); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            void top_loop(tags::Bool<true>) {}

            void top_loop(tags::Bool<false>) {
                try { level_loop(tags::Ind<1>()); }
                catch (const ExitLoop &ex) {throw ex;}
            }

        public:
            Ordered(size_t n) :
                    Base<nind>(nterm(n)), m_n(n) {}

        protected:
            void throwing_loop() override {
                top_loop(tags::Bool<nind == 0>());
            }
        };
    }

    /**
     * "run time number of dimensions"
     */
    namespace rtnd {
        using inds_t = std::vector<size_t>;

        struct Base {
        protected:
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t m_value;
            /**
             * integer index of the current iteration (body call)
             */
            size_t m_iiter;
        public:
            /**
             * number of dimensions: length of the index array
             */
            const size_t m_nind;
            /**
             * length of the index iterator
             */
            const size_t m_niter;

            const inds_t &value() const { return m_value; }

            const size_t &iiter() const { return m_iiter; }

            size_t ind_sum() const {
                return std::accumulate(m_value.cbegin(), m_value.cend(), 0ul);
            }

            void set_null_value() {
                m_value.assign(m_nind, ~0ul);
            }

            Base(size_t nind, size_t nterm) : m_value(nind, 0), m_nind(nind), m_niter(nterm) {}

            /**
             * function to be called each time a new set of indices is formed
             */
            virtual void body() = 0;

        protected:
            /**
             * function with defines the looping logic, calls body, and increments the iteration counter
             */
            virtual void throwing_loop() = 0;

        public:
            /**
             * exposed wrapper function which catches any instance of the ExitLoop exception thrown, and gracefully
             * handles the zero-dimensional edge case
             */
            void loop() {
                if (m_nind) {
                    m_iiter = ~0ul;
                    try {
                        throwing_loop();
                        // if loop completes without early termination, increment to final value
                        ++m_iiter;
                        DEBUG_ASSERT_EQ(m_iiter, m_niter, "loop completed after incorrect number of iterations");
                    }
                    catch (const ExitLoop &) {}
                } else m_iiter = 0ul;
                // leave value "undefined"
                set_null_value();
            }
        };

        struct Unrestricted : Base {
            const inds_t m_shape;
        private:
            static size_t nterm(const inds_t &shape) {
                if (shape.empty()) return 0;
                size_t n = 1;
                for (const auto &i: shape) n *= i;
                return n;
            }

            void level_loop(size_t ilevel) {
                const auto &iind = ilevel - 1;
                auto &ind = m_value[iind];
                const auto &extent = m_shape[iind];
                try {
                    if (ilevel < m_nind)
                        for (ind = 0ul; ind < extent; ++ind) level_loop(ilevel + 1);
                    else {
                        for (ind = 0ul; ind < extent; ++ind) {
                            ++m_iiter;
                            body();
                        }
                    }
                }
                catch (const ExitLoop& ex){throw ex;}
            }
        public:
            Unrestricted(const inds_t &shape) : Base(shape.size(), nterm(shape)), m_shape(shape) {}

            Unrestricted(const size_t &nind, const size_t &extent) : Unrestricted(std::vector<size_t>(nind, extent)) {}

        protected:
            void throwing_loop() override {
                level_loop(1);
            }
        };

        template<bool strict = true, bool ascending = true>
        struct Ordered : Base {
        private:
            const size_t m_n;

            static size_t nterm(size_t n, size_t r) {
                if (!r) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + r) - 1, r);
            }

            void level_loop(size_t ilevel) {
                const size_t iind = ascending ? (m_nind - ilevel) : (ilevel - 1);
                const size_t iind_unrestrict = ascending ? m_nind - 1 : 0;
                auto &ind = m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                try {
                    if (ilevel < m_nind)
                        for (ind = 0ul; ind < extent; ++ind) level_loop(ilevel + 1);
                    else {
                        for (ind = 0ul; ind < extent; ++ind) {
                            ++m_iiter;
                            body();
                        }
                    }
                }
                catch (const ExitLoop& ex){throw ex;}
            }

        public:
            Ordered(const size_t &n, const size_t &r) :
                    Base(r, nterm(n, r)), m_n(n) {}

        protected:
            void throwing_loop() override {
                level_loop(1);
            }
        };
    }
}


#endif //M7_FOREACHVIRTUAL_H
