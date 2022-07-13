//
// Created by Robert J. Anderson on 16/06/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/linalg/Sparse.h"

TEST(SparseMatrix, RealMultiplication){
    typedef double T;
    const uint_t n=3;
    v_t<T> v(n, 0);
    sparse::dynamic::Matrix<T> matrix;
    for(uint_t i=0; i<n; ++i) {
        v[i] = i+1;
        for(uint_t j=0; j<n; ++j) {
            matrix.add(j, i, i*3.0-j*2.0+1);
        }
    }
    v_t<T> res(n, 0);
    matrix.multiply(v, res);
    ASSERT_EQ(res[0], 30);
    ASSERT_EQ(res[1], 18);
    ASSERT_EQ(res[2], 6);
}

TEST(SparseMatrix, ComplexMultiplication){
    typedef std::complex<double> T;
    const uint_t n=3;
    v_t<T> v(n, 0);
    sparse::dynamic::Matrix<T> matrix;
    for(uint_t i=0; i<n; ++i) {
        v[i] = T{double(i+1), 2};
        for(uint_t j=0; j<n; ++j) {
            matrix.add(j, i, T{i*3.0-j*2.0+1, double(i+j)});
        }
    }
    v_t<T> res(n, 0);
    matrix.multiply(v, res);
    ASSERT_EQ(res[0].real(), 24);
    ASSERT_EQ(res[0].imag(), 32);
    ASSERT_EQ(res[1].real(), 6);
    ASSERT_EQ(res[1].imag(), 26);
    ASSERT_EQ(res[2].real(), -12);
    ASSERT_EQ(res[2].imag(), 20);
}