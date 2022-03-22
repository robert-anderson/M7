//
// Created by rja on 10/01/2022.
//

#include "Examples.h"
#include "M7_lib/hash/Hashing.h"
#include "M7_lib/parallel/MPIAssert.h"

sparse::Matrix<double>
sparse_matrix_examples::rect_double(const size_t &nrow, const size_t &ncol, const size_t &nnonzero_per_row) {
    REQUIRE_LE(nnonzero_per_row, ncol,
               "specified number of non-zero elements per row exceeds column count");
    const size_t v_hi = 7;
    sparse::Matrix<double> out;
    out.resize(nrow);

    defs::inds icols;
    std::vector<double> values;

    for (size_t irow = 0ul; irow < nrow; ++irow) {
        utils::convert(hashing::unique_in_range(irow, nnonzero_per_row, 0, ncol, true), icols);
        utils::convert(hashing::in_range(irow, nnonzero_per_row, 1, v_hi, true), values);
        out.insert(irow, icols, values);
    }

    return out;
}