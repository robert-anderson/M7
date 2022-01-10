//
// Created by anderson on 1/10/22.
//

#ifndef M7_SPARSE_EXAMPLES_H
#define M7_SPARSE_EXAMPLES_H

#include "src/core/parallel/MPIAssert.h"

namespace sparse_matrix_examples {

    sparse::Matrix<double> rect_double(const size_t &nrow, const size_t &ncol, const size_t &nnonzero_per_row) {
        REQUIRE_LE(nnonzero_per_row, ncol,
                        "specified number of non-zero elements per row exceeds column count");
        const size_t v_hi = 4;
        sparse::Matrix<double> out;
        out.resize(nrow);

        defs::inds icols;
        std::vector<double> values;

        for (size_t irow = 0ul; irow < nrow; ++irow) {
            utils::convert(hashing::unique_in_range(irow, nnonzero_per_row, 0, ncol, true), icols);
            utils::convert(hashing::in_range(irow, nnonzero_per_row, 1, v_hi, true), values);
            out.add(irow, icols, values);
        }

        return out;
    }


}

#endif //M7_SPARSE_EXAMPLES_H
