//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSEMATRIX_H
#define M7_SPARSEMATRIX_H

#include <vector>
#include <forward_list>
#include <src/core/parallel/MPIAssert.h>
#include "src/core/linalg/Matrix.h"
#include "src/core/io/Logging.h"

namespace sparse {

    class Network {
        bool m_resized_by_add = false;

    protected:
        std::vector<std::forward_list<size_t>> m_rows_icols;

    public:
        void resize(const size_t nrow);

        void expand(const size_t delta_nrow);

        size_t nrow() const;

        void add(const size_t &irow, const size_t &icol);

        bool empty();

        const std::forward_list<size_t>& operator[](const size_t& irow);
    };

    template<typename T>
    class Matrix : public Network {
        std::vector<std::forward_list<T>> m_rows_values;

    public:
        void add(const size_t &irow, const size_t &icol, const T& v) {
            Network::add(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            m_rows_values[irow].push_front(v);
        }

        void multiply(const std::vector<T> &in, std::vector<T> &out) const {
            // each row of the out vector receives a contribution from every
            // corresponding entry in the sparse matrix row
            // out_i = sum_j H_ij in_j
            for (size_t irow = 0; irow < m_rows_icols.size(); ++irow) {
                auto icol_it = m_rows_icols[irow].cbegin();
                auto value_it = m_rows_values[irow].cbegin();
                while (icol_it!=m_rows_icols[irow].cend()){
                    DEBUG_ASSERT_LT(*icol_it, in.size(), "sparse matrix column OOB");
                    DEBUG_ASSERT_FALSE(value_it==m_rows_values[irow].cend(),
                                       "values list incongruent with column indices list");
                    out[irow] += *value_it * in[*icol_it];
                    ++icol_it;
                    ++value_it;
                }
            }
        }

        std::pair<const std::forward_list<size_t>&, const std::forward_list<T>&> operator[](const size_t& irow) {
            ASSERT(irow<nrow());
            return {m_rows_icols[irow], m_rows_values[irow]};
        }
    };
};

#endif //M7_SPARSEMATRIX_H
