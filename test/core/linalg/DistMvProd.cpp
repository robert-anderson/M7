//
// Created by anderson on 1/7/22.
//

#include <src/core/hash/Hashing.h>
#include "src/core/linalg/DistMvProd.h"
#include "test/core/sparse/Examples.h"
#include "gtest/gtest.h"

TEST(DistMvProd, SparseRealSym) {
    typedef double T;
    const size_t nrow = 123;
    const size_t ncol = 19;
    const size_t nnonzero_per_row = 5;
    auto global_mat = sparse_matrix_examples::rect_double(nrow, ncol, nnonzero_per_row);
    const size_t count_local = mpi::evenly_shared_count(nrow);
    const size_t displ_local = mpi::evenly_shared_displ(nrow);
    auto local_mat = global_mat.get_row_subset(count_local, displ_local);

    std::vector<T> in, out, out_chk;
    if (mpi::i_am_root()) utils::convert(hashing::in_range(0, nrow, 0, 12, false), in);

    dist_mv_prod::Sparse<T> prod(local_mat);
    prod.parallel_multiply(in, out);

    if (mpi::i_am_root()) {
        global_mat.multiply(in, out_chk);
        ASSERT_EQ(out, out_chk);
    }
}