//
// Created by anderson on 2/9/22.
//

#ifndef M7_FOREACHVIRTUAL_H
#define M7_FOREACHVIRTUAL_H

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
             * work space into which each set of indices is inserted: should not be directly accessed
             */
            inds_t<nind> m_inds;
        public:
            /**
             * getter for the current index value so that body implementations make immutable reference to the indices
             */
            const inds_t<nind> &inds() const { return m_inds; }
            /**
             * length of the index iterator
             */
            const size_t m_nterm;

            /**
             * @return
             *  sum over current indices
             */
            size_t sum() const {
                return std::accumulate(m_inds.cbegin(), m_inds.cend(), 0);
            }

            Base(size_t nterm) : m_nterm(nterm) {}

            /**
             * function to be called each time a new set of indices is formed
             */
            virtual void body() = 0;

            /**
             * function defining the algorithm by which sets of indices are generated
             */
            virtual void loop() = 0;

        };

        template<size_t nind>
        struct Unrestricted : Base<nind> {
            using Base<nind>::m_inds;
            using Base<nind>::inds;
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
             * loop through the values of the ilevel element of m_inds between 0 and the extent of that dimension
             * @tparam ilevel
             *  current level (position in the m_inds array)
             */
            template<size_t ilevel>
            void level_loop(tags::Ind<ilevel>) {
                constexpr size_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_inds[iind];
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
                auto &ind = m_inds[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    try { body(); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            /**
             * in the edge case that the nind is 0, do nothing
             */
            void top_loop(tags::Bool<true>) {}

            /**
             * if nind is nonzero, start at the first index
             */
            void top_loop(tags::Bool<false>) {
                try { level_loop(tags::Ind<1>()); }
                catch (const ExitLoop &) {}
            }

        public:
            Unrestricted(const inds_t<nind> &shape) : Base<nind>(nterm(shape)), m_shape(shape) {}

            Unrestricted(const size_t &extent) : Unrestricted(array_utils::filled<size_t, nind>(extent)) {}

            void loop() override {
                top_loop(tags::Bool<nind == 0>());
            }
        };

        template<size_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_inds;
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
                auto &ind = Base<nind>::m_inds[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    try { level_loop(tags::Ind<ilevel + 1>()); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            void level_loop(tags::Ind<nind>) {
                constexpr size_t iind = ascending ? 0 : nind - 1;
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_inds[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    try { body(); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            void top_loop(tags::Bool<true>) {}

            void top_loop(tags::Bool<false>) {
                try { level_loop(tags::Ind<1>()); }
                catch (const ExitLoop &) {}
            }

        public:
            Ordered(size_t n) :
                    Base<nind>(nterm(n)), m_n(n) {}

            void loop() override {
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
            /**
             * work space into which each set of indices is inserted: should not be directly accessed
             */
            inds_t m_inds;
            /**
             * number of dimensions: length of the index array
             */
            const size_t m_nind;
            /**
             * length of the index iterator
             */
            const size_t m_nterm;

            const inds_t &inds() const { return m_inds; }

            size_t ind_sum() const {
                return std::accumulate(m_inds.cbegin(), m_inds.cend(), 0ul);
            }

            Base(size_t nind, size_t nterm) : m_inds(nind, 0), m_nind(nind), m_nterm(nterm) {}

            virtual void body() = 0;

            virtual void loop() = 0;

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
                auto &ind = m_inds[iind];
                const auto &extent = m_shape[iind];
                try {
                    if (ilevel < m_nind)
                        for (ind = 0ul; ind < extent; ++ind) level_loop(ilevel + 1);
                    else
                        for (ind = 0ul; ind < extent; ++ind) body();
                }
                catch (const ExitLoop& ex){throw ex;}
            }
        public:
            Unrestricted(const inds_t &shape) : Base(shape.size(), nterm(shape)), m_shape(shape) {}

            Unrestricted(const size_t &nind, const size_t &extent) : Unrestricted(std::vector<size_t>(nind, extent)) {}

            void loop() override {
                if (m_nind == 0) return;
                try {level_loop(1);}
                catch (const ExitLoop&) {}
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
                auto &ind = m_inds[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1] + !strict;
                try {
                    if (ilevel < m_nind)
                        for (ind = 0ul; ind < extent; ++ind) level_loop(ilevel + 1);
                    else
                        for (ind = 0ul; ind < extent; ++ind) body();
                }
                catch (const ExitLoop& ex){throw ex;}
            }

        public:
            Ordered(const size_t &n, const size_t &r) :
                    Base(r, nterm(n, r)), m_n(n) {}

            void loop() override {
                if (m_nind == 0) return;
                try {level_loop(1);}
                catch (const ExitLoop&){}
            }
        };
    }
}


#endif //M7_FOREACHVIRTUAL_H
