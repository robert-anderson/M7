//
// Created by Robert J. Anderson on 1/10/22.
//

#ifndef M7_SPARSE_EXAMPLES_H
#define M7_SPARSE_EXAMPLES_H

#include "M7_lib/linalg/Sparse.h"

namespace sparse_matrix_examples {

    sparse::Matrix<double> rect_double(const size_t &nrow, const size_t &ncol, const size_t &nnonzero_per_row);
    sparse::Matrix<std::complex<double>> rect_double_complex(const size_t &nrow, const size_t &ncol, const size_t &nnonzero_per_row);

}

#endif //M7_SPARSE_EXAMPLES_H
