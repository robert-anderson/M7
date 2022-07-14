//
// Created by Robert J. Anderson on 1/7/22.
//

#include <M7_lib/util/Hash.h>
#include <M7_lib/linalg/DistMvProd.h>
#include <test_core/sparse/Examples.h>
#include "gtest/gtest.h"

TEST(DistMvProd, SparseRealSym) {
    typedef double T;
    const uint_t nrow = 123;
    const uint_t ncol = 19;
    const uint_t nnonzero_per_row = 5;
    auto global_mat = sparse_matrix_examples::rect_double(nrow, ncol, nnonzero_per_row);
    const uint_t count_local = mpi::evenly_shared_count(nrow);
    const uint_t displ_local = mpi::evenly_shared_displ(nrow);
    auto local_mat = global_mat.row_subset(count_local, displ_local);
    ASSERT_EQ(global_mat.nrow(), nrow);
    ASSERT_EQ(mpi::all_sum(local_mat.nrow()), nrow);

    v_t<T> in, out, out_chk;
    if (mpi::i_am_root()) convert::vector(hash::in_range(0, nrow, 0, 12, false), in);

    dist_mv_prod::Sparse<T> prod(local_mat);
    prod.parallel_multiply(in, out, false);

    if (mpi::i_am_root()) {
        global_mat.multiply(in, out_chk);
        ASSERT_EQ(out, out_chk);
    }
}
