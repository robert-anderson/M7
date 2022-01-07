//
// Created by anderson on 1/7/22.
//

#include <src/core/hash/Hashing.h>
#include "src/core/linalg/DistMvProd.h"
#include "gtest/gtest.h"

TEST(DistMvProd, SparseRealSym) {
    typedef double T;
    sparse::Matrix<T> global_mat;
    const size_t nrow = 40;
    const size_t ncol = 24;
    const size_t nnonzero_per_row = 6;
    const size_t v_hi=10;
    for (size_t irow=0ul; irow<nrow; ++irow) {
        auto icols = hashing::in_range(irow, nnonzero_per_row, 0, ncol, true);
        auto vs = hashing::in_range(irow, nnonzero_per_row, 0, v_hi, false);
        global_mat.add(irow, icols, vs);
    }
}