//
// Created by Robert J. Anderson on 10/01/2022.
//

#include "Examples.h"
#include "M7_lib/util/Hash.h"
#include "M7_lib/parallel/MPIAssert.h"

sparse::dynamic::Matrix<double>
sparse_matrix_examples::rect_double(const uint_t &nrow, const uint_t &ncol, const uint_t &nnonzero_per_row) {
    REQUIRE_LE(nnonzero_per_row, ncol,
               "specified number of non-zero elements per row exceeds column count");
    const uint_t v_hi = 7;
    sparse::dynamic::Matrix<double> out;
    out.resize(nrow);

    uintv_t icols;
    v_t<double> values;

    for (uint_t irow = 0ul; irow < nrow; ++irow) {
        convert::vector(hash::unique_in_range<>(irow, nnonzero_per_row, 0, ncol, true), icols);
        convert::vector(hash::in_range(irow, nnonzero_per_row, 1, v_hi, true), values);
        out.insert(irow, icols, values);
    }

    return out;
}

sparse::dynamic::Matrix<std::complex<double>>
sparse_matrix_examples::rect_double_complex(const uint_t &nrow, const uint_t &ncol, const uint_t &nnonzero_per_row) {
    auto real = rect_double(nrow, ncol, nnonzero_per_row);
    sparse::dynamic::Matrix<std::complex<double>> out;
    out.resize(nrow);
    for (uint_t irow=0ul; irow<nrow; ++irow) {
        const auto real_row = real[irow];
        for (uint_t ielem=0ul; ielem<real_row.first.size(); ++ielem){
            double imag = hash::in_range({irow, ielem}, 0, 4);
            out.add(irow, {real_row.first[ielem], {real_row.second[ielem], imag}});
        }
    }
    return out;
}
