//
// Created by Robert J. Anderson on 14/06/2020.
//

#ifndef M7_SPARSE_H
#define M7_SPARSE_H

#include <vector>
#include <forward_list>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/io/Logging.h>
#include <set>
#include "M7_lib/util/Convert.h"

namespace sparse {

    struct Element {
        /**
         * sparsely-connected index
         */
        uint_t m_i;
        virtual str_t to_string() const {
            return convert::to_string(m_i);
        }
    };

    template<typename T>
    struct MatrixElement : Element {
        /**
         * value associated with the sparsely-connected element
         */
        T m_v;

        str_t to_string() const override {
            return log::format("({} -> {})", m_i, convert::to_string(m_v));
        }
    };

    namespace dynamic {
        struct Base {
        protected:
            const uint_t m_init_row_nentry;
            uint_t m_nentry = 0ul;
            uint_t m_max_col_ind = 0ul;

            Base(uint_t init_row_nentry) : m_init_row_nentry(init_row_nentry) {}

            virtual str_t row_to_string(uint_t irow) const = 0;

        public:
            virtual str_t to_string() const = 0;
        };


        template<typename T>
        struct Generic : Base {
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
        private:
            typedef v_t<v_t<T>> rows_t;
            rows_t m_rows;
        public:

            Generic(uint_t nrow = 0ul, uint_t init_row_nentry = 6ul) :
                    Base(init_row_nentry), m_rows(nrow, v_t<T>(init_row_nentry)) {
                for (auto &row: m_rows) row.clear();
            }

            /*
             * ctor for a subset block of rows
             */
            Generic(const Generic &other, uint_t count, uint_t displ) :
                    Generic(count, other.m_init_row_nentry) {
                REQUIRE_LE(displ, nrow(), "row offset OOB");
                REQUIRE_LE(displ + count, nrow(), "row offset+count OOB");
                auto begin = m_rows.cbegin();
                std::advance(begin, displ);
                auto end = begin;
                std::advance(end, count);
                m_rows = rows_t(begin, end);
                // data is now copied, now update metadata
                for (const auto &row: m_rows) {
                    if (row.empty()) continue;
                    m_nentry += row.size();
                    auto max_element = std::max_element(row.cbegin(), row.cend());
                    m_max_col_ind = std::max(m_max_col_ind, *max_element);
                }
            }

            Generic &operator=(const Generic &other) {
                if (&other == this) return *this;
                m_rows = other.m_rows;
                m_max_col_ind = other.m_max_col_ind;
                m_nentry = other.m_nentry;
                return *this;
            }

            Generic(const Generic &other, bool symmetrize) : Generic(other.nrow(), other.m_init_row_nentry) {
                if (!symmetrize) {
                    *this = other;
                    return;
                }
                REQUIRE_LT(other.m_max_col_ind, nrow(), "too many columns for this to be a symmetric matrix");
                for (uint_t irow = 0ul; irow < nrow(); ++irow) {
                    const auto &row = m_rows[irow];
                    for (uint_t ientry = 0ul; ientry < row.size(); ++ientry) {
                        const auto &entry = row[ientry];
                        add(irow, entry);
                        if (entry != irow) insert(entry, irow);
                    }
                }
            }

            Generic(const Generic &other) : Generic(other, false) {}

            /**
             * reallocate the m_rows vector
             * @param nrow
             *  new length of m_rows
             */
            void resize(uint nrow) {
                const auto nrow_old = m_rows.size();
                if (nrow <= nrow_old) return;
                m_rows.resize(nrow, v_t<T>(m_init_row_nentry));
                for (uint_t irow = nrow_old; irow < nrow; ++irow) m_rows[irow].clear();
            }

            uint_t nrow() const {
                return m_rows.size();
            }

            const v_t<T> &operator[](uint_t irow) const {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return m_rows[irow];
            }

        protected:
            /*
             * no modification of the rows outside of this class definition
             */
            v_t<T> &operator[](uint_t irow) {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return m_rows[irow];
            }

        public:
            /**
             * append the element to the given row
             * @param irow
             *  row index onto which the element is appended
             * @param elem
             *  element to append
             */
            void add(uint_t irow, const T &elem) {
                if (irow >= nrow()) resize(irow + 1);
                auto base_elem = static_cast<const Element &>(elem);
                (*this)[irow].emplace_back(elem);
                const auto icol = base_elem.m_i;
                if (icol > m_max_col_ind) m_max_col_ind = icol;
                ++m_nentry;
                //return m_rows[irow].size()-1;
            }

            /**
             * lookup the element index in the row first. if it exists, overwrite it, else simply delegate to add
             * @param irow
             *  row index
             * @param elem
             *  element which is either added, or clobbers existing element
             */
            void insert(uint_t irow, const T &elem) {
                auto &row = (*this)[irow];
                const auto icol = static_cast<const Element &>(elem).m_i;
                auto pred = [&icol](const T &v) {
                    return static_cast<const Element &>(v).m_i == icol;
                };
                auto it = std::find_if(row.begin(), row.end(), pred);
                if (it == row.end()) {
                    add(irow, elem);
                } else {
                    *it = elem;
                }
            }

            bool empty() const {
                return !m_nentry;
            }

            bool empty(uint_t irow) const {
                return (*this)[irow].empty();
            }

            str_t to_string() const override {
                str_t out;
                for (uint_t irow = 0ul; irow < nrow(); ++irow) {
                    if (m_rows[irow].empty()) continue;
                    out += std::to_string(irow) + ": " + row_to_string(irow) + "\n";
                }
                return out;
            }

        protected:
            str_t row_to_string(uint_t irow) const override {
                strv_t out;
                const auto &row = m_rows[irow];
                for (auto &entry: row)
                    out.push_back(static_cast<const Element &>(entry).to_string());
                return convert::to_string(out);
            }

        };

        typedef Generic<Element> Network;
        template<typename T>
        using Matrix = Generic<MatrixElement<T>>;
    }

    /**
     * number of rows is variable, and the vectors storing the contents of each can be resized. Ideal for building up
     * sparse structures from an empty state
     */
    namespace dynamic {

        class Network {
            bool m_resized_by_add = false;
            uint_t m_nentry = 0ul;
            uint_t m_max_icol = 0ul;

        protected:
            v_t<uintv_t> m_rows_icols;

        public:

            virtual void resize(uint_t nrow);

            uint_t nrow() const;

            uint_t nentry() const;

            uint_t nentry(uint_t irow) const;

            uint_t max_column_index() const;

            uint_t add(uint_t irow, uint_t icol);

            uint_t insert(uint_t irow, uint_t icol);

            void add(uint_t irow, const uintv_t &icols);

            void insert(uint_t irow, const uintv_t &icols);

            bool empty() const;

            bool empty(uint_t irow) const;

            const uintv_t &operator[](uint_t irow) const;

            virtual str_t row_to_string(uint_t irow) const;

            str_t to_string() const;

            Network get_symmetrized() const;

            Network get_row_subset(uint_t count, uint_t displ) const;

        protected:
            void get_row_subset(Network &subnet, uint_t count, uint_t displ) const;
        };

        template<typename T>
        class Matrix : public sparse::dynamic::Network {
            v_t<v_t<T>> m_rows_values;

        public:

            void resize(uint_t nrow) override {
                Network::resize(nrow);
                if (nrow > m_rows_values.size()) m_rows_values.resize(nrow);
            }

            uint_t add(uint_t irow, uint_t icol, const T &v) {
                auto i = Network::add(irow, icol);
                if (irow >= m_rows_values.size()) resize(irow + 1);
                m_rows_values[irow].push_back(v);
                return i;
            }

            uint_t insert(uint_t irow, uint_t icol, const T &v) {
                auto i = Network::insert(irow, icol);
                if (irow >= m_rows_values.size()) resize(irow + 1);
                auto& row = m_rows_values[irow];
                if (i < row.size()) row[i] = v;
                else row.push_back(v);
                return i;
            }

            void add(uint_t irow, const uintv_t &icols, const v_t<T> &vs) {
                REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
                for (uint_t i = 0ul; i < icols.size(); ++i) add(irow, icols[i], vs[i]);
            }

            void insert(uint_t irow, const uintv_t &icols, const v_t<T> &vs) {
                REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
                for (uint_t i = 0ul; i < icols.size(); ++i) insert(irow, icols[i], vs[i]);
            }

            void add(uint_t irow, const std::pair<uint_t, T> &pair) {
                add(irow, pair.first, pair.second);
            }

            void insert(uint_t irow, const std::pair<uint_t, T> &pair) {
                insert(irow, pair.first, pair.second);
            }

            void multiply(const T *v, T *mv, uint_t mv_size) const {
                // each row of the mv vector receives a contribution from every
                // corresponding entry v the sparse matrix row
                // (Mv)_i = sum_j M_ij v_j
                memset(static_cast<void*>(mv), 0, sizeof(T) * mv_size);
                for (uint_t irow = 0; irow < m_rows_icols.size(); ++irow) {
                    auto icol_it = m_rows_icols[irow].cbegin();
                    auto value_it = m_rows_values[irow].cbegin();
                    for (; icol_it != m_rows_icols[irow].cend(); (++icol_it, ++value_it)) {
                        DEBUG_ASSERT_FALSE(value_it == m_rows_values[irow].cend(),
                                           "values list incongruent with column indices list");
                        mv[irow] += *value_it * v[*icol_it];
                    }
                }
            }

            void multiply(const v_t<T> &v, v_t<T> &mv) const {
                if (nrow() > mv.size()) mv.resize(nrow());
                multiply(v.data(), mv.data(), mv.size());
            }

            std::pair<const uintv_t &, const v_t<T> &> operator[](uint_t irow) const {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return {m_rows_icols[irow], m_rows_values[irow]};
            }

            str_t row_to_string(uint_t irow) const override {
                strv_t out;
                for (uint_t ientry = 0ul; ientry < m_rows_icols[irow].size(); ++ientry) {
                    uint_t icol = m_rows_icols[irow][ientry];
                    auto v = m_rows_values[irow][ientry];
                    out.push_back("(" + convert::to_string(icol) + " -> " + convert::to_string(v) + ")");
                }
                return convert::to_string(out);
            }

            Matrix<T> get_symmetrized(bool conj) const {
                Matrix<T> sym_mat;
                sym_mat.resize(nrow());
                REQUIRE_LT(max_column_index(), nrow(), "too many columns for this to be a symmetric matrix");
                for (uint_t irow=0ul; irow<nrow(); ++irow) {
                    const auto& icols = m_rows_icols[irow];
                    const auto& values = m_rows_values[irow];
                    for (uint_t iicol=0ul; iicol < icols.size(); ++iicol) {
                        const auto& icol = icols[iicol];
                        const auto& value = values[iicol];
                        sym_mat.insert(irow, icol, value);
                        if (icol != irow) sym_mat.insert(icol, irow, conj ? arith::conj(value): value);
                    }
                }
                return sym_mat;
            }

            Matrix<T> get_row_subset(uint_t count, uint_t displ) const {
                Matrix<T> submat;
                Network::get_row_subset(submat, count, displ);
                auto begin = m_rows_values.cbegin()+uintv_t::difference_type(displ);
                auto end = begin + uintv_t::difference_type(count);
                REQUIRE_GE(std::distance(end, m_rows_values.cend()), 0, "end iterator OOB");
                submat.m_rows_values = v_t<v_t<T>>(begin, end);
                return submat;
            }
        };
    }

    /**
     * number of entries is constant, therefore all data can be stored contiguously in memory. This should confer
     * performance benefits
     */
    namespace fixed {

        class Base {

            static uintv_t make_counts(const dynamic::Network& src);

            static uintv_t make_displs(const uintv_t& counts);

        public:
            const uint_t m_nrow;
            const uint_t m_max_nentry;
        protected:
            const uintv_t m_displs;
            const uint_t m_nentry;
            Base(const uintv_t& counts);

            Base(const dynamic::Network& src);
        };

        class Network : public Base {
            const uintv_t m_entries;

            static uintv_t make_entries(const dynamic::Network& src);

            uintv_t::const_iterator citer(uint_t irow) const;

        public:
            explicit Network(const dynamic::Network& src);

            uintv_t::const_iterator cbegin(uint_t irow) const;

            uintv_t::const_iterator cend(uint_t irow) const;

            uint_t nentry(uint_t irow) const;
        };

        template<typename T>
        class Matrix : public Base {
        public:
            typedef std::pair<uint_t, T> pair_t;
            typedef v_t<pair_t> entries_t;
        private:
            const entries_t m_entries;

            entries_t make_entries(const dynamic::Matrix<T>& src) {
                v_t<pair_t> out;
                out.reserve(src.nrow());
                for (uint_t irow=0ul; irow<src.nrow(); ++irow) {
                    const auto row = convert::zip<uint_t, int>(src[irow]);
                    out.insert(out.cend(), row.cbegin(), row.cend());
                }
                return out;
            }

            typename entries_t::const_iterator citer(uint_t irow) const {
                auto it = m_entries.cbegin();
                std::advance(it, m_displs[irow]);
                return it;
            }

        public:
            Matrix(const dynamic::Matrix<T>& src): Base(src), m_entries(make_entries(src)){
                DEBUG_ASSERT_EQ(m_entries.size(), m_nentry, "incorrect number of entries");
            }

            typename entries_t::const_iterator cbegin(uint_t irow) const {
                DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
                return citer(irow);
            }

            typename entries_t::const_iterator cend(uint_t irow) const {
                DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
                return citer(irow+1);
            }

            uint_t nentry(uint_t irow) const {
                return std::distance(cbegin(irow), cend(irow));
            }

            const pair_t& get(uint_t irow, uint_t ientry) const {
                DEBUG_ASSERT_LT(ientry, nentry(irow), "entry index OOB");
                auto it = cbegin(irow);
                std::advance(it, ientry);
                return *it;
            }
        };
    }

    /**
     * more space-efficient representations of the inverse mappings of sparse structures than are available from the
     * dense::Matrix
     */
    namespace inverse {

        class Network {
            std::set<std::pair<uint_t, uint_t>> m_conns;
        public:
            Network(const dynamic::Network& src);

            Network(const fixed::Network& src);

            bool lookup(uint_t i, uint_t j) const;
        };

        template<typename T>
        class Matrix {
            std::map<std::pair<uint_t, uint_t>, T> m_conns;
        public:
            Matrix(const dynamic::Matrix<T>& src) {
                for (uint_t irow=0ul; irow<src.nrow(); ++irow) {
                    const auto& inds = src[irow].first;
                    const auto& values = src[irow].second;
                    for (uint_t i=0ul; i<inds.size(); ++i) m_conns.insert({irow, inds[i]}, values[i]);
                }
            }

            Matrix(const fixed::Matrix<T>& src) {
                for (auto irow=0ul; irow<src.m_nrow; ++irow)
                    for (auto it = src.cbegin(irow); it != src.cend(irow); ++it)
                        m_conns.insert({{irow, it->first}, it->second});
            }

            /**
             * first element in returned pair is true only if the entry was found in the map
             */
            std::pair<bool, T> lookup(uint_t i, uint_t j) const {
                auto it = m_conns.find({i, j});
                return {it!=m_conns.cend(), it!=m_conns.cend() ? it->second : T{}};
            }
        };
    }
}

#endif //M7_SPARSE_H
