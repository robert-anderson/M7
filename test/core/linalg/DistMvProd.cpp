//
// Created by anderson on 1/7/22.
//

#include <src/core/hash/Hashing.h>
#include "src/core/linalg/DistMvProd.h"
#include "gtest/gtest.h"

TEST(DistMvProd, SparseRealSym) {
    typedef double T;
    const size_t nrow = 50;
    const size_t ncol = 43;
    const size_t nnonzero_per_row = 6;
    const size_t v_hi=4;
    sparse::Matrix<T> global_mat;
    global_mat.resize(nrow);
    sparse::Matrix<T> local_mat;
    const size_t nrow_local = mpi::evenly_shared_count(nrow);
    const size_t displ_local = mpi::evenly_shared_displ(nrow);
    local_mat.resize(nrow_local);

    defs::inds icols;
    std::vector<T> values;

    if (mpi::i_am_root()) {
        for (size_t irow = 0ul; irow < nrow; ++irow) {
            utils::convert(hashing::unique_in_range(irow, nnonzero_per_row, 0, ncol, true), icols);
            utils::convert(hashing::in_range(irow, nnonzero_per_row, 1, v_hi, true), values);
            global_mat.add(irow, icols, values);
        }
    }

    for (size_t irow = 0ul; irow < nrow_local; ++irow) {
        auto irow_global = irow + displ_local;
        utils::convert(hashing::unique_in_range(irow_global, nnonzero_per_row, 0, ncol, true), icols);
        utils::convert(hashing::in_range(irow_global, nnonzero_per_row, 1, v_hi, true), values);
        local_mat.add(irow, icols, values);
    }

    std::vector<T> in, out, out_chk;
    if (mpi::i_am_root()) utils::convert(hashing::in_range(0, nrow, 0, 12, false), in);

    dist_mv_prod::Sparse<T> prod(local_mat);
    prod.parallel_multiply(in, out);

    if (mpi::i_am_root()) {
        global_mat.multiply(in, out_chk);
        ASSERT_EQ(out, out_chk);
    }
}