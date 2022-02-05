//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSE_H
#define M7_SPARSE_H

#include <vector>
#include <forward_list>
#include <src/core/parallel/MPIAssert.h>
#include "src/core/io/Logging.h"

namespace sparse {

    class Network {
        bool m_resized_by_add = false;
        size_t m_nentry = 0ul;
        size_t m_max_icol = 0ul;

    protected:
        std::vector<defs::inds> m_rows_icols;

    public:

        virtual void resize(const size_t &nrow);

        size_t nrow() const;

        size_t nentry() const;

        size_t nentry(const size_t& irow) const;

        size_t max_column_index() const;

        size_t add(const size_t &irow, const size_t &icol);

        size_t insert(const size_t &irow, const size_t &icol);

        void add(const size_t &irow, const defs::inds &icols);

        void insert(const size_t &irow, const defs::inds &icols);

        bool empty() const;

        bool empty(const size_t& irow) const;

        const defs::inds &operator[](const size_t &irow) const;

        virtual std::vector<std::string> row_to_strings(size_t irow) const;

        std::string to_string() const;

        Network get_symmetrized() const;

        Network get_row_subset(size_t count, size_t displ) const;

    protected:
        void get_row_subset(Network &subnet, size_t count, size_t displ) const;
    };

    template<typename T>
    class Matrix : public Network {
        std::vector<std::vector<T>> m_rows_values;

    public:

        void resize(const size_t &nrow) override {
            Network::resize(nrow);
            if (nrow > m_rows_values.size()) m_rows_values.resize(nrow);
        }

        size_t add(const size_t &irow, const size_t &icol, const T &v) {
            auto i = Network::add(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            m_rows_values[irow].push_back(v);
            return i;
        }

        size_t insert(const size_t &irow, const size_t &icol, const T &v) {
            auto i = Network::insert(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            auto& row = m_rows_values[irow];
            if (i < row.size()) row[i] = v;
            else row.push_back(v);
            return i;
        }

        void add(const size_t &irow, const defs::inds &icols, const std::vector<T> &vs) {
            REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
            for (size_t i = 0ul; i < icols.size(); ++i) add(irow, icols[i], vs[i]);
        }

        void insert(const size_t &irow, const defs::inds &icols, const std::vector<T> &vs) {
            REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
            for (size_t i = 0ul; i < icols.size(); ++i) insert(irow, icols[i], vs[i]);
        }

        void add(const size_t &irow, const std::pair<size_t, T> &pair) {
            add(irow, pair.first, pair.second);
        }

        void insert(const size_t &irow, const std::pair<size_t, T> &pair) {
            insert(irow, pair.first, pair.second);
        }

        void multiply(const T *v, T *mv, size_t mv_size) const {
            // each row of the mv vector receives a contribution from every
            // corresponding entry v the sparse matrix row
            // (Mv)_i = sum_j M_ij v_j
            std::memset(static_cast<void*>(mv), 0, sizeof(T) * mv_size);
            for (size_t irow = 0; irow < m_rows_icols.size(); ++irow) {
                auto icol_it = m_rows_icols[irow].cbegin();
                auto value_it = m_rows_values[irow].cbegin();
                for (; icol_it != m_rows_icols[irow].cend(); (++icol_it, ++value_it)) {
                    DEBUG_ASSERT_FALSE(value_it == m_rows_values[irow].cend(),
                                       "values list incongruent with column indices list");
                    mv[irow] += *value_it * v[*icol_it];
                }
            }
        }

        void multiply(const std::vector<T> &v, std::vector<T> &mv) const {
            if (nrow() > mv.size()) mv.resize(nrow());
            multiply(v.data(), mv.data(), mv.size());
        }

        std::pair<const defs::inds &, const std::vector<T> &> operator[](const size_t &irow) const {
            DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
            return {m_rows_icols[irow], m_rows_values[irow]};
        }

        std::vector<std::string> row_to_strings(size_t irow) const override {
            std::vector<std::string> out;
            for (size_t ientry = 0ul; ientry < m_rows_icols[irow].size(); ++ientry) {
                size_t icol = m_rows_icols[irow][ientry];
                auto v = m_rows_values[irow][ientry];
                out.push_back("(" + utils::to_string(icol) + " -> " + utils::to_string(v) + ")");
            }
            return out;
        }

        Matrix<T> get_symmetrized(bool conj) const {
            Matrix<T> sym_mat;
            sym_mat.resize(nrow());
            REQUIRE_LT(max_column_index(), nrow(), "too many columns for this to be a symmetric matrix");
            for (size_t irow=0ul; irow<nrow(); ++irow) {
                const auto& icols = m_rows_icols[irow];
                const auto& values = m_rows_values[irow];
                for (size_t iicol=0ul; iicol < icols.size(); ++iicol) {
                    const auto& icol = icols[iicol];
                    const auto& value = values[iicol];
                    sym_mat.insert(irow, icol, value);
                    if (icol != irow) sym_mat.insert(icol, irow, value);
                }
            }
            return sym_mat;
        }

        Matrix<T> get_row_subset(size_t count, size_t displ) const {
            Matrix<T> submat;
            Network::get_row_subset(submat, count, displ);
            auto begin = m_rows_values.cbegin()+defs::inds::difference_type(displ);
            auto end = begin + defs::inds::difference_type(count);
            REQUIRE_GE(std::distance(end, m_rows_values.cend()), 0, "end iterator OOB");
            submat.m_rows_values = std::vector<std::vector<T>>(begin, end);
            return submat;
        }
    };
}

#endif //M7_SPARSE_H
