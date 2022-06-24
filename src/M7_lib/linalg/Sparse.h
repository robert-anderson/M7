//
// Created by Robert J. Anderson on 14/06/2020.
//

#ifndef M7_SPARSE_H
#define M7_SPARSE_H

#include <vector>
#include <forward_list>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/io/Logging.h>
#include "M7_lib/util/Convert.h"

namespace sparse {

    class Network {
        bool m_resized_by_add = false;
        uint_t m_nentry = 0ul;
        uint_t m_max_icol = 0ul;

    protected:
        std::vector<uintv_t> m_rows_icols;

    public:

        virtual void resize(const uint_t &nrow);

        uint_t nrow() const;

        uint_t nentry() const;

        uint_t nentry(const uint_t& irow) const;

        uint_t max_column_index() const;

        uint_t add(const uint_t &irow, const uint_t &icol);

        uint_t insert(const uint_t &irow, const uint_t &icol);

        void add(const uint_t &irow, const uintv_t &icols);

        void insert(const uint_t &irow, const uintv_t &icols);

        bool empty() const;

        bool empty(const uint_t& irow) const;

        const uintv_t &operator[](const uint_t &irow) const;

        virtual std::string row_to_string(uint_t irow) const;

        std::string to_string() const;

        Network get_symmetrized() const;

        Network get_row_subset(uint_t count, uint_t displ) const;

    protected:
        void get_row_subset(Network &subnet, uint_t count, uint_t displ) const;
    };

    template<typename T>
    class Matrix : public Network {
        std::vector<std::vector<T>> m_rows_values;

    public:

        void resize(const uint_t &nrow) override {
            Network::resize(nrow);
            if (nrow > m_rows_values.size()) m_rows_values.resize(nrow);
        }

        uint_t add(const uint_t &irow, const uint_t &icol, const T &v) {
            auto i = Network::add(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            m_rows_values[irow].push_back(v);
            return i;
        }

        uint_t insert(const uint_t &irow, const uint_t &icol, const T &v) {
            auto i = Network::insert(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            auto& row = m_rows_values[irow];
            if (i < row.size()) row[i] = v;
            else row.push_back(v);
            return i;
        }

        void add(const uint_t &irow, const uintv_t &icols, const std::vector<T> &vs) {
            REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
            for (uint_t i = 0ul; i < icols.size(); ++i) add(irow, icols[i], vs[i]);
        }

        void insert(const uint_t &irow, const uintv_t &icols, const std::vector<T> &vs) {
            REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
            for (uint_t i = 0ul; i < icols.size(); ++i) insert(irow, icols[i], vs[i]);
        }

        void add(const uint_t &irow, const std::pair<uint_t, T> &pair) {
            add(irow, pair.first, pair.second);
        }

        void insert(const uint_t &irow, const std::pair<uint_t, T> &pair) {
            insert(irow, pair.first, pair.second);
        }

        void multiply(const T *v, T *mv, uint_t mv_size) const {
            // each row of the mv vector receives a contribution from every
            // corresponding entry v the sparse matrix row
            // (Mv)_i = sum_j M_ij v_j
            std::memset(static_cast<void*>(mv), 0, sizeof(T) * mv_size);
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

        void multiply(const std::vector<T> &v, std::vector<T> &mv) const {
            if (nrow() > mv.size()) mv.resize(nrow());
            multiply(v.data(), mv.data(), mv.size());
        }

        std::pair<const uintv_t &, const std::vector<T> &> operator[](const uint_t &irow) const {
            DEBUG_ASSERT_LT(irow, nrow(), "row index OOB");
            return {m_rows_icols[irow], m_rows_values[irow]};
        }

        std::string row_to_string(uint_t irow) const override {
            std::vector<std::string> out;
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
            submat.m_rows_values = std::vector<std::vector<T>>(begin, end);
            return submat;
        }
    };
}

#endif //M7_SPARSE_H
