//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSE_H
#define M7_SPARSE_H

#include <vector>
#include <forward_list>
#include <src/core/parallel/MPIAssert.h>
#include "src/core/linalg/Matrix.h"
#include "src/core/io/Logging.h"

namespace sparse {

    class Network {
        bool m_resized_by_add = false;

    protected:
        std::vector<defs::inds> m_rows_icols;

    public:
        virtual void resize(const size_t& nrow);

        size_t nrow() const;

        void add(const size_t &irow, const size_t &icol);

        bool empty();

        const defs::inds& operator[](const size_t& irow);
    };

    template<typename T>
    class Matrix : public Network {
        std::vector<std::vector<T>> m_rows_values;

    public:

        void resize(const size_t &nrow) override {
            Network::resize(nrow);
            if (nrow > m_rows_values.size()) m_rows_values.resize(nrow);
        }

        void add(const size_t &irow, const size_t &icol, const T& v) {
            Network::add(irow, icol);
            if (irow >= m_rows_values.size()) resize(irow + 1);
            m_rows_values[irow].push_back(v);
        }

        void add(const size_t &irow, const std::pair<size_t, T>& pair){
            add(irow, pair.first, pair.second);
        }

        void multiply(const std::vector<T> &in, std::vector<T> &out) const {
            // each row of the out vector receives a contribution from every
            // corresponding entry in the sparse matrix row
            // out_i = sum_j H_ij in_j
            for (size_t irow = 0; irow < m_rows_icols.size(); ++irow) {
                auto icol_it = m_rows_icols[irow].cbegin();
                auto value_it = m_rows_values[irow].cbegin();
                for (; icol_it!=m_rows_icols[irow].cend(); (++icol_it, ++value_it)){
                    DEBUG_ASSERT_LT(*icol_it, in.size(), "sparse matrix column OOB");
                    DEBUG_ASSERT_FALSE(value_it==m_rows_values[irow].cend(),
                                       "values list incongruent with column indices list");
                    out[irow] += *value_it * in[*icol_it];
                }
            }
        }

        std::pair<const defs::inds&, const std::vector<T>&> operator[](const size_t& irow) const {
            ASSERT(irow<nrow());
            return {m_rows_icols[irow], m_rows_values[irow]};
        }
    };
};

#endif //M7_SPARSE_H
