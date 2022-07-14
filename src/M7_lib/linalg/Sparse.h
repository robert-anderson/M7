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

    /**
     * In a sparse structure, rows are made up of many elements.
     * each row corresponds to a single basis vector of a matrix, or node of a network
     * A network element carries no more than an integer identifying the connected node index.
     */
    struct Element {
        /**
         * sparsely-connected index
         */
        uint_t m_i;
        Element(uint_t i=~0ul): m_i(i){}

        virtual ~Element() = default;

        virtual str_t to_string() const {
            return convert::to_string(m_i);
        }

        operator uint_t () const {
            return m_i;
        }

        template<typename T>
        static const Element& cast(const T& elem) {
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            return elem;
        }

        template<typename T>
        static Element& cast(T& elem) {
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            return elem;
        }
    };

    template<typename T>
    struct MatrixElement : Element {
        /**
         * value associated with the sparsely-connected element
         */
        T m_v;

        MatrixElement(uint_t i=~0ul, T v={}): Element(i), m_v(v){}

        str_t to_string() const override {
            return log::format("({} -> {})", m_i, convert::to_string(m_v));
        }
    };

    /**
     * number of rows is variable, and the vectors storing the contents of each can be resized. Ideal for building up
     * sparse structures from an empty state
     */
    namespace dynamic {
        struct Base {
        protected:
            const uint_t m_init_row_nentry;
            uint_t m_nentry = 0ul;
            uint_t m_max_col_ind = 0ul;

            Base(uint_t init_row_nentry) : m_init_row_nentry(init_row_nentry) {}
            virtual ~Base() = default;

            virtual str_t row_to_string(uint_t irow) const = 0;

        public:
            bool empty() const {
                return !m_nentry;
            }

            uint_t nentry() const {
                return m_nentry;
            }

            virtual uint_t nentry(uint_t irow) const = 0;

            virtual uint_t nrow() const = 0;

            bool empty(uint_t irow) const {
                return !nentry(irow);
            }

            uint_t max_col_ind() const {
                return m_max_col_ind;
            }

            virtual str_t to_string() const = 0;

        };


        template<typename T>
        struct Generic : Base {
            typedef T elem_t;
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
        protected:
            typedef v_t<T> row_t;
            typedef v_t<row_t> rows_t;
            rows_t m_rows;

            template<typename fn_t>
            void set_symmetrized(const Generic& src, fn_t fn) {
                functor::assert_prototype<T(T, uint_t)>(fn);
                resize(src.nrow());
                REQUIRE_LT(src.max_col_ind(), src.nrow(), "too many columns for this to be a symmetric matrix");
                for (uint_t irow = 0ul; irow < nrow(); ++irow) {
                    const auto& row = src.get(irow);
                    for (auto& elem : row) {
                        uint_t icol = Element::cast(elem);
                        insert(irow, elem);
                        if (icol != irow) insert(icol, fn(elem, irow));
                    }
                }
            }

        public:

            Generic(uint_t nrow = 0ul, uint_t init_row_nentry = 6ul) :
                    Base(init_row_nentry), m_rows(nrow, row_t(init_row_nentry, T())) {
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
                /*
                 * predicate to compare column indices so that m_max_col_ind can be updated properly
                 */
                auto pred = [](const Element& e1, const Element& e2){
                    return e1.m_i > e2.m_i;
                };
                for (const auto &row: m_rows) {
                    if (row.empty()) continue;
                    m_nentry += row.size();
                    auto max_element = std::max_element(row.cbegin(), row.cend(), pred);
                    m_max_col_ind = std::max(m_max_col_ind, Element::cast(*max_element).m_i);
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
                auto fn = [](T elem, uint_t irow) {
                    auto out = elem;
                    Element::cast(out).m_i = irow;
                    return out;
                };
                set_symmetrized(other, fn);
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
                m_rows.resize(nrow, row_t(m_init_row_nentry));
                for (uint_t irow = nrow_old; irow < nrow; ++irow) m_rows[irow].clear();
            }

        protected:
            /*
             * no modification of the rows outside of this class definition and descendents
             */
            row_t & get(uint_t irow) {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return m_rows[irow];
            }

            const row_t & get(uint_t irow) const {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return m_rows[irow];
            }

        public:

            const row_t &operator[](uint_t irow) const {
                return get(irow);
            }

            /**
             * append the element to the given row
             * @param irow
             *  row index onto which the element is appended
             * @param elem
             *  element to append
             */
            void add(uint_t irow, const T &elem) {
                if (irow >= nrow()) resize(irow + 1);
                get(irow).emplace_back(elem);
                const auto icol = Element::cast(elem).m_i;
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
                auto &row = get(irow);
                const auto icol = Element::cast(elem).m_i;
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

            str_t to_string() const override {
                str_t out;
                for (uint_t irow = 0ul; irow < nrow(); ++irow) {
                    if (m_rows[irow].empty()) continue;
                    out += std::to_string(irow) + ": " + row_to_string(irow) + "\n";
                }
                return out;
            }

            uint_t nentry(uint_t irow) const override {
                return (*this)[irow].size();
            }

            uint_t nrow() const override {
                return m_rows.size();
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
        struct Matrix : Generic<MatrixElement<T>> {
            typedef MatrixElement<T> elem_t;
            using Generic<elem_t>::nrow;
            using Generic<elem_t>::get;
            using Generic<elem_t>::set_symmetrized;
            using Generic<elem_t>::insert;
            Matrix(uint_t nrow = 0ul, uint_t init_row_nentry = 6ul): Generic<elem_t>(nrow, init_row_nentry){}

            Matrix(const Matrix &other, uint_t count, uint_t displ) : Generic<elem_t>(other, count, displ){}

            Matrix(const Matrix &other, bool conj) : Generic<elem_t>(other.nrow(), other.m_init_row_nentry) {
                auto fn = [&conj](elem_t elem, uint_t irow) {
                    return elem_t(irow, conj ? arith::conj(elem.m_v) : elem.m_v);
                };
                set_symmetrized(other, fn);
            }

            Matrix(const Matrix &other) : Generic<elem_t>(other, false) {}

            void insert(uint_t irow, const uintv_t& icols, const v_t<T>& values) {
                REQUIRE_EQ(icols.size(), values.size(), "number of column indices and values inserted must match");
                for (uint_t i=0ul; i<icols.size(); ++i) insert(irow, {icols[i], values[i]});
            }

            Matrix symmetrized(bool conj) {
                return {*this, conj};
            }

            Matrix row_subset(uint_t count, uint_t displ) {
                return {*this, count, displ};
            }

            void multiply(const T *v, T *mv, uint_t mv_size) const {
                // each row of the mv vector receives a contribution from every
                // corresponding entry v the sparse matrix row
                // (Mv)_i = sum_j M_ij v_j
                memset(static_cast<void*>(mv), 0, sizeof(T) * mv_size);
                for (uint_t irow = 0; irow < nrow(); ++irow) {
                    const auto& row = get(irow);
                    for (const elem_t& elem : row) mv[irow] += elem.m_v * v[elem.m_i];
                }
            }

            void multiply(const v_t<T> &v, v_t<T> &mv) const {
                if (nrow() > mv.size()) mv.resize(nrow());
                multiply(v.data(), mv.data(), mv.size());
            }
        };
    }

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
    }

    /**
     * more space-efficient representations of the inverse mappings of sparse structures than are available from the
     * dense::Matrix
     */
    namespace inverse {

        template<typename T>
        struct Generic {
            typedef T elem_t;
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            typedef std::pair<uint_t, uint_t> key_t;
            typedef std::pair<key_t, T> pair_t;
        private:
            const T m_not_found_entry;
            std::map<key_t, T> m_conns;
        public:
            Generic(const dynamic::Generic<T>& src) {
                for (uint_t irow=0ul; irow<src.nrow(); ++irow) {
                    for (auto& it: src[irow]){
                        const auto elem = static_cast<const Element&>(it);
                        m_conns.insert({{irow, elem.m_i}, it});
                    }
                }
            }

            Generic(const fixed::Generic<T>& src) {
                for (uint_t irow=0ul; irow<src.m_nrow; ++irow) {
                    for (auto it = src.cbegin(irow); it != src.cend(irow); ++it) {
                        const auto elem = static_cast<const Element&>(*it);
                        m_conns.insert({{irow, elem.m_i}, *it});
                    }
                }
            }

            const T& entry(uint_t i, uint_t j) const {
                auto it = m_conns.find({i, j});
                return it==m_conns.cend() ? m_not_found_entry : it->second;
            }

            bool exists(uint_t i, uint_t j) const {
                return m_conns.find({i, j}) != m_conns.cend();
            }
        };

        typedef Generic<Element> Network;

        template<typename T>
        struct Matrix : Generic<MatrixElement<T>> {
            typedef MatrixElement<T> elem_t;
            Matrix(const dynamic::Matrix<T>& src): Generic<elem_t>(src){}
            Matrix(const fixed::Matrix<T>& src): Generic<elem_t>(src){}

            T get(uint_t i, uint_t j) const {
                return Generic<elem_t>::entry(i, j).m_v;
            }
        };
    }

#if 0
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
#endif
}

#endif //M7_SPARSE_H
