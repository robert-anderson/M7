//
// Created by rja on 14/07/22.
//

#ifndef M7_SPARSE_DYNAMIC_H
#define M7_SPARSE_DYNAMIC_H

#include "Sparse.h"

namespace sparse {
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
            typedef std::list<T> row_t;
            typedef v_t<row_t> rows_t;
            rows_t m_rows;

            template<typename fn_t>
            void set_symmetrized(const Generic& src, fn_t fn) {
                functor::assert_prototype<T(T, uint_t)>(fn);
                resize(src.nrow());
                REQUIRE_LT(src.max_col_ind(), src.nrow(), "too many columns for this to be a symmetric matrix");
                for (uint_t irow = 0ul; irow < nrow(); ++irow) {
                    const auto& row = src.get(irow);
                    for (auto& elem: row) {
                        uint_t icol = Element::cast(elem);
                        insert(irow, elem);
                        if (icol != irow) insert(icol, fn(elem, irow));
                    }
                }
            }

        public:

            Generic(uint_t nrow = 0ul, uint_t init_row_nentry = 6ul) :
                    Base(init_row_nentry), m_rows(nrow, row_t(init_row_nentry, T())) {
                for (auto& row: m_rows) row.clear();
            }

            /*
             * ctor for a subset block of rows
             */
            Generic(const Generic& other, uint_t count, uint_t displ) :
                    Generic(count, other.m_init_row_nentry) {
                REQUIRE_LE(displ, other.nrow(), "row offset OOB");
                REQUIRE_LE(displ + count, other.nrow(), "row offset+count OOB");
                auto begin = other.m_rows.cbegin();
                std::advance(begin, displ);
                auto end = begin;
                std::advance(end, count);
                m_rows = rows_t(begin, end);
                // data is now copied, now update metadata
                /*
                 * predicate to compare column indices so that m_max_col_ind can be updated properly
                 */
                auto pred = [](const Element& e1, const Element& e2) {
                    return e1.m_i > e2.m_i;
                };
                for (const auto& row: m_rows) {
                    if (row.empty()) continue;
                    m_nentry += row.size();
                    auto max_element = std::max_element(row.cbegin(), row.cend(), pred);
                    m_max_col_ind = std::max(m_max_col_ind, Element::cast(*max_element).m_i);
                }
            }

            Generic& operator=(const Generic& other) {
                if (&other == this) return *this;
                m_rows = other.m_rows;
                m_max_col_ind = other.m_max_col_ind;
                m_nentry = other.m_nentry;
                return *this;
            }

            Generic(const Generic& other, bool symmetrize) : Generic(other.nrow(), other.m_init_row_nentry) {
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

            Generic(const Generic& other) : Generic(other, false) {}

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
            row_t& get(uint_t irow) {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return m_rows[irow];
            }

            const row_t& get(uint_t irow) const {
                DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
                return m_rows[irow];
            }

        private:
            T* find(uint_t irow, uint_t icol) {
                auto& row = get(irow);
                if (row.empty()) return nullptr;
                auto pred = [&icol](const T& v) {
                    return static_cast<const Element&>(v).m_i == icol;
                };
                auto it = std::find_if(row.begin(), row.end(), pred);
                return it == row.end() ? nullptr : &(*it);
            }

        public:

            const row_t& operator[](uint_t irow) const {
                return get(irow);
            }

            /**
             * append the element to the given row
             * @param irow
             *  row index onto which the element is appended
             * @param elem
             *  element to append
             */
            void add(uint_t irow, const T& elem) {
                if (irow >= nrow()) resize(irow + 1);
                const auto icol = Element::cast(elem).m_i;
                DEBUG_ASSERT_FALSE(find(irow, icol), "multiple occurrences of the same connected index in row");
                get(irow).emplace_back(elem);
                if (icol > m_max_col_ind) m_max_col_ind = icol;
                ++m_nentry;
            }

            /**
             * lookup the element index in the row first. if it exists, overwrite it, else simply delegate to add
             * @param irow
             *  row index
             * @param elem
             *  element which is either added, or clobbers existing element
             */
            void insert(uint_t irow, const T& elem) {
                if (irow >= nrow()) {
                    add(irow, elem);
                    return;
                }
                auto ptr = find(irow, Element::cast(elem).m_i);
                if (ptr) *ptr = elem;
                else add(irow, elem);
            }

            bool operator==(const Generic<T>& other) const {
                if (nrow()!=other.nrow()) return false;
                if (m_nentry!=other.m_nentry) return false;
                if (m_max_col_ind!=other.m_max_col_ind) return false;
                return m_rows==other.m_rows;
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
                const auto& row = m_rows[irow];
                for (auto& entry: row)
                    out.push_back(static_cast<const Element&>(entry).to_string());
                return string::join(out, ", ");
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

            Matrix(uint_t nrow = 0ul, uint_t init_row_nentry = 6ul) : Generic<elem_t>(nrow, init_row_nentry) {}

            Matrix(const Matrix& other, uint_t count, uint_t displ) : Generic<elem_t>(other, count, displ) {}

            Matrix(const Matrix& other, bool conj) : Generic<elem_t>(other.nrow(), other.m_init_row_nentry) {
                auto fn = [&conj](elem_t elem, uint_t irow) {
                    return elem_t(irow, conj ? arith::conj(elem.m_v) : elem.m_v);
                };
                set_symmetrized(other, fn);
            }

            Matrix(const Matrix& other) : Generic<elem_t>(other, false) {}

            void insert(uint_t irow, const uintv_t& icols, const v_t<T>& values) {
                REQUIRE_EQ(icols.size(), values.size(), "number of column indices and values inserted must match");
                for (uint_t i = 0ul; i < icols.size(); ++i) insert(irow, {icols[i], values[i]});
            }

            Matrix symmetrized(bool conj) {
                return {*this, conj};
            }

            Matrix row_subset(uint_t count, uint_t displ) {
                return {*this, count, displ};
            }

            void multiply(const T* v, T* mv, uint_t mv_size) const {
                // each row of the mv vector receives a contribution from every
                // corresponding entry v the sparse matrix row
                // (Mv)_i = sum_j M_ij v_j
                memset(static_cast<void*>(mv), 0, sizeof(T) * mv_size);
                for (uint_t irow = 0; irow < nrow(); ++irow) {
                    const auto& row = get(irow);
                    for (const elem_t& elem: row) mv[irow] += elem.m_v * v[elem.m_i];
                }
            }

            void multiply(const v_t<T>& v, v_t<T>& mv) const {
                if (nrow() > mv.size()) mv.resize(nrow());
                multiply(v.data(), mv.data(), mv.size());
            }
        };
    } // dynamic
} // sparse

#endif //M7_SPARSE_DYNAMIC_H
