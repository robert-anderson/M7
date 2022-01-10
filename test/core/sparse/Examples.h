//
// Created by anderson on 1/10/22.
//

#ifndef M7_SPARSE_EXAMPLES_H
#define M7_SPARSE_EXAMPLES_H

#include "src/core/parallel/MPIAssert.h"

namespace sparse_matrix_examples {

    template<typename T>
    struct DistPair {
        /**
         * stored on this rank only
         */
        sparse::Matrix<T> m_local;
        /**
         * the whole matrix across all ranks (as though m_local were gathered row-wise)
         */
        sparse::Matrix<T> m_global;

        DistPair(const size_t &nrow, const defs::inds& irows,
                 const defs::inds& icols, const std::vector<T>& values) {
            REQUIRE_EQ(irows.size(), icols.size(), "row and column index vectors must have equal size");
            REQUIRE_EQ(irows.size(), values.size(), "row index and value vectors must have equal size");
        }

        DistPair<double> dist_double(const size_t &nrow, const size_t &ncol,
                                     const size_t &nnonzero_per_row, bool global_root_only = true) {
            DEBUG_ASSERT_LE(nonzero_per_row, ncol,
                            "specified number of non-zero elements per row exceeds column count");
            bool gets_global = mpi::i_am_root() || !global_root_only;
            const size_t v_hi = 4;
            DistPair<double> out;
            if (gets_global) out.m_global.resize(nrow);
            const size_t nrow_local = mpi::evenly_shared_count(nrow);
            const size_t displ_local = mpi::evenly_shared_displ(nrow);
            out.m_local.resize(nrow_local);

            defs::inds icols;
            std::vector <T> values;

            if (gets_global) {
                for (size_t irow = 0ul; irow < nrow; ++irow) {
                    utils::convert(hashing::unique_in_range(irow, nnonzero_per_row, 0, ncol, true), icols);
                    utils::convert(hashing::in_range(irow, nnonzero_per_row, 1, v_hi, true), values);
                    out.m_global.add(irow, icols, values);
                }
            }

            for (size_t irow = 0ul; irow < nrow_local; ++irow) {
                auto irow_global = irow + displ_local;
                utils::convert(hashing::unique_in_range(irow_global, nnonzero_per_row, 0, ncol, true), icols);
                utils::convert(hashing::in_range(irow_global, nnonzero_per_row, 1, v_hi, true), values);
                out.m_local.add(irow, icols, values);
            }
            return out;
        }


    }

#endif //M7_SPARSE_EXAMPLES_H
