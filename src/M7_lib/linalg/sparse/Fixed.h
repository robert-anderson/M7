//
// Created by rja on 14/07/22.
//

#ifndef M7_SPARSE_FIXED_H
#define M7_SPARSE_FIXED_H

#include "Dynamic.h"

namespace sparse {
    /**
     * number of entries is constant, therefore all data can be stored contiguously in memory. This should confer
     * performance benefits
     */
    namespace fixed {

        struct Base {
            const uint_t m_nrow;
            const uint_t m_max_nentry;
        protected:
            const uintv_t m_displs;
            const uint_t m_nentry;

            Base(const uintv_t& counts);

            Base(const dynamic::Base& src);

            static uintv_t make_counts(const dynamic::Base& src);

            static uintv_t make_displs(const uintv_t& counts);
        };


        template<typename T>
        struct Generic : Base {
            typedef T elem_t;
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
        private:
            typedef v_t<T> entries_t;
            const v_t<T> m_entries;

            static entries_t make_entries(const dynamic::Generic<T>& src) {
                entries_t out;
                out.reserve(src.nrow());
                for (uint_t irow = 0ul; irow < src.nrow(); ++irow)
                    out.insert(out.cend(), src[irow].cbegin(), src[irow].cend());
                return out;
            }

        public:

            Generic(const dynamic::Generic<T>& src) : Base(src), m_entries(make_entries(src)) {
                DEBUG_ASSERT_EQ(m_entries.size(), m_nentry, "incorrect number of entries");
            }

            typename entries_t::const_iterator cbegin(uint_t irow) const {
                if (irow >= m_nrow) return m_entries.cend();
                auto it = m_entries.cbegin();
                std::advance(it, m_displs[irow]);
                return it;
            }

            typename entries_t::const_iterator cend(uint_t irow) const {
                return cbegin(irow + 1);
            }

            uint_t nentry(uint_t irow) const {
                DEBUG_ASSERT_LT(irow, m_nrow, "Row index OOB");
                return m_displs[irow + 1] - m_displs[irow];
            }

            const T& get(uint_t irow, uint_t ientry) const {
                DEBUG_ASSERT_LT(ientry, nentry(irow), "entry index OOB");
                auto it = cbegin(irow);
                std::advance(it, ientry);
                return *it;
            }
        };

        typedef Generic<Element> Network;
        template<typename T>
        using Matrix = Generic<MatrixElement<T>>;
    } // fixed
} // sparse

#endif //M7_SPARSE_FIXED_H
