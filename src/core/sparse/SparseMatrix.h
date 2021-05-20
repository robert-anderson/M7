//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSEMATRIX_H
#define M7_SPARSEMATRIX_H

#include <vector>
#include <forward_list>
#include "src/core/linalg/Matrix.h"
#include "src/core/io/Logging.h"
#include "omp.h"

namespace sparse {
    template<typename T>
    struct Entry {
        size_t icol;
        T element;
    };

    template<typename T>
    struct ParallelVector {
        const size_t m_local_size;
        const defs::inds m_counts, m_displs;
        const size_t m_global_size;
        std::vector<T> m_local;
        std::vector<T> m_global;

        ParallelVector(size_t local_size) :
        m_local_size(local_size),
        m_counts(mpi::all_gathered(m_local_size)), m_displs(mpi::counts_to_displs_consec(m_counts)),
        m_global_size(m_displs.back()+m_counts.back()), m_local(m_local_size, 0){
            m_global.resize(m_global_size, 0);
        }

        void all_gatherv() {
            mpi::all_gatherv(m_local.data(), m_local_size, m_global.data(), m_counts, m_displs);
        }
    };

    template<typename T>
    class Matrix {
        bool m_resized_by_add = false;
        std::vector<std::forward_list<Entry<T>>> m_rows;

    public:
        void resize(const size_t nrow) {
            ASSERT(nrow >= m_rows.size());
            m_rows.resize(nrow);
        }

        void expand(const size_t delta_nrow) {
            resize(m_rows.size() + delta_nrow);
        }

        size_t nrow() const {
            return m_rows.size();
        }

        void add(const size_t &irow, const size_t &icol, const T& v) {
            if (irow >= m_rows.size()) {
                if (!m_resized_by_add) {
                    log::warn("Resizing SparseMatrix by adding a row (this entails reallocation which is inefficient)");
                    log::warn("Call resize before filling if number of rows is known in advance");
                    m_resized_by_add = true;
                }
                resize(irow + 1);
            }
            m_rows[irow].push_front(Entry<T>{icol, v});
        }

        bool empty() { return m_rows.empty(); }

        void multiply(const std::vector<T> &in, std::vector<T> &out) const {
            // each row of the out vector receives a contribution from every
            // corresponding entry in the sparse matrix row
            // out_i = sum_j H_ij in_j
            for (size_t irow = 0; irow < m_rows.size(); ++irow) {
                for (auto entry = m_rows[irow].begin(); entry != m_rows[irow].end(); entry++) {
                    auto &jrow = entry->icol;
                    ASSERT(jrow < in.size())
                    auto &hij = entry->element;
                    ASSERT(hij == hij)
                    ASSERT(in[jrow] == in[jrow])
                    out[irow] += hij * in[jrow];
                    ASSERT(out[irow] == out[irow])
                }
            }
        }

        void multiply(ParallelVector<T>& v){
            multiply(v.m_global, v.m_local);
            v.all_gatherv();
        }

        const std::forward_list<Entry<T>> row(const size_t& irow) {
            ASSERT(irow<nrow());
            return m_rows[irow];
        }
    };


};

#endif //M7_SPARSEMATRIX_H
