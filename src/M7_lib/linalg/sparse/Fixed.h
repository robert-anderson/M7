//
// Created by rja on 14/07/22.
//

#ifndef M7_SPARSE_FIXED_H
#define M7_SPARSE_FIXED_H

#include "M7_lib/hdf5/File.h"
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
            virtual ~Base() = default;

            static uintv_t make_counts(const dynamic::Base& src);

            static uintv_t make_displs(const uintv_t& counts);

            virtual str_t row_to_string(uint_t irow) const = 0;

        public:
            virtual str_t to_string() const = 0;
        };


        template<typename T>
        struct Generic : Base {
            typedef T elem_t;
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            /**
             * true if all rows entries are in ascending order of column index
             */
            const bool m_ordered;
        private:
            typedef v_t<T> entries_t;
            const v_t<T> m_entries;

            static entries_t make_entries(const dynamic::Generic<T>& src, bool ordered) {
                entries_t out;
                out.reserve(src.nrow());
                for (uint_t irow = 0ul; irow < src.nrow(); ++irow) {
                    out.insert(out.end(), src[irow].cbegin(), src[irow].cend());
                    if (ordered) {
                        /*
                         * operate only on the block just inserted
                         */
                        auto end = out.end();
                        auto begin = end;
                        std::advance(begin, -src[irow].size());
                        std::sort(begin, end);
                    }
                }
                return out;
            }

        public:

            Generic(const dynamic::Generic<T>& src, bool ordered=true) :
                Base(src), m_ordered(ordered), m_entries(make_entries(src, ordered)) {
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

            bool operator==(const Generic<T>& other) const {
                if (m_nrow!=other.m_nrow) return false;
                if (m_nentry!=other.m_nentry) return false;
                if (m_max_nentry!=other.m_max_nentry) return false;
                return m_entries==other.m_entries;
            }

            str_t to_string() const override {
                str_t out;
                for (uint_t irow = 0ul; irow < m_nrow; ++irow) {
                    if (!nentry(irow)) continue;
                    out += std::to_string(irow) + ": " + row_to_string(irow) + "\n";
                }
                return out;
            }

        private:

            // TODO: update with new HDF5 saves
//            void write_column_inds(hdf5::NodeWriter& parent, const v_t<T>& entries, uint_t irank=0ul) const {
//                uintv_t inds;
//                inds.reserve(entries.size());
//                for (auto& entry : entries) inds.push_back(static_cast<const Element&>(entry).m_i);
//                parent.write_data("col_indices", inds, irank);
//            }
//
//            void write_values(hdf5::NodeWriter& parent, const v_t<Element>& entries, uint_t irank=0ul) const {
//                // Element base class has no values
//            }
//
//            template<typename U>
//            void write_values(hdf5::NodeWriter& parent, const v_t<MatrixElement<U>>& entries, uint_t irank=0ul) const {
//                v_t<U> values;
//                values.reserve(entries.size());
//                for (auto& entry : entries) values.push_back(entry.m_v);
//                parent.write_data("values", values, irank);
//            }
//
//        public:
//            void save(hdf5::NodeWriter& parent, uint_t irank=0ul) const {
//                v_t<uint_t> nentries;
//                nentries.reserve(m_nrow);
//                v_t<uint_t> icols;
//                icols.reserve(m_nentry);
//                v_t<T> values;
//                values.reserve(m_nentry);
//                for (uint_t irow=0ul; irow<m_nrow; ++irow) {
//                    nentries.push_back(nentry(irow));
//                }
//                parent.write_data("row_sizes", nentries, irank);
//                write_column_inds(parent, m_entries, irank);
//                write_values(parent, m_entries, irank);
//            }

        protected:
            str_t row_to_string(uint_t irow) const override {
                strv_t out;
                for (auto it = cbegin(irow); it != cend(irow); ++it)
                    out.push_back(static_cast<const Element&>(*it).to_string());
                return string::join(out, ", ");
            }
        };

        typedef Generic<Element> Network;
        template<typename T>
        using Matrix = Generic<MatrixElement<T>>;
    } // fixed
} // sparse

#endif //M7_SPARSE_FIXED_H
