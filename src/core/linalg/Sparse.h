//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSE_H
#define M7_SPARSE_H

#include <vector>
#include <forward_list>
#include <src/core/parallel/MPIAssert.h>
#include "src/core/linalg/Dense.h"
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

        size_t max_column_index() const;

        void add(const size_t &irow, const size_t &icol);

        void checked_add(const size_t &irow, const size_t &icol);

        void add(const size_t &irow, const defs::inds &icols);

        bool empty();

        const defs::inds &operator[](const size_t &irow);

        virtual std::vector<std::string> row_to_strings(size_t irow) const;

        std::string to_string() const;

        Network get_symmetrized() const;
    };

    template<typename T>
    class Matrix : public Network {
        std::vector<std::vector<T>> m_rows_values;

    public:

        void resize(const size_t &nrow) override {
            Network::resize(nrow);
            if (nrow > m_rows_values.size()) m_rows_values.resize(nrow);
        }

        void add(const size_t &irow, const size_t &icol, const T &v) {
            Network::add(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            m_rows_values[irow].push_back(v);
        }

        void add(const size_t &irow, const defs::inds &icols, const std::vector<T> &vs) {
            REQUIRE_EQ(icols.size(), vs.size(), "must have same number of column indices and values");
            for (size_t i = 0ul; i < icols.size(); ++i) add(irow, icols[i], vs[i]);
        }

        void add(const size_t &irow, const std::pair<size_t, T> &pair) {
            add(irow, pair.first, pair.second);
        }

        void multiply(const T *in, T *out) const {
            // each row of the out vector receives a contribution from every
            // corresponding entry in the sparse matrix row
            // out_i = sum_j H_ij in_j
            for (size_t irow = 0; irow < m_rows_icols.size(); ++irow) {
                auto icol_it = m_rows_icols[irow].cbegin();
                auto value_it = m_rows_values[irow].cbegin();
                for (; icol_it != m_rows_icols[irow].cend(); (++icol_it, ++value_it)) {
                    DEBUG_ASSERT_FALSE(value_it == m_rows_values[irow].cend(),
                                       "values list incongruent with column indices list");
                    out[irow] += *value_it * in[*icol_it];
                }
            }
        }

        void multiply(const std::vector<T> &in, std::vector<T> &out) const {
            if (nrow() > out.size()) out.resize(nrow());
            multiply(in.data(), out.data());
        }

        std::pair<const defs::inds &, const std::vector<T> &> operator[](const size_t &irow) const {
            ASSERT(irow < nrow());
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

        dense::Matrix<T> to_dense() const {
            dense::Matrix<T> mat(nrow());
            for (size_t irow = 0ul; irow < nrow(); ++irow) {
                for (size_t ientry = 0ul; ientry < m_rows_icols[irow].size(); ++ientry) {
                    size_t icol = m_rows_icols[irow][ientry];
                    auto v = m_rows_values[irow][ientry];
                    mat(irow, icol) = v;
                }
            }
            return mat;
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
                    sym_mat.add(irow, icol, value);
                    if (icol != irow) sym_mat.checked_add(icol, irow, value);
                }
            }
            return sym_mat;
        }
    };
}

#endif //M7_SPARSE_H
