//
// Created by rja on 16/06/2020.
//

#include <gtest/gtest.h>
#include "src/core/linalg/Matrix.h"

TEST(Matrix, SquareRealMultiplication){
    typedef double T;
    const size_t n=3;
    std::vector<T> v(n, 0);
    Matrix<T> matrix(n);
    for(size_t i=0; i<n; ++i) {
        v[i] = i+1;
        for(size_t j=0; j<n; ++j) {
            matrix(i, j) = i*3.0-j*2.0+1;
        }
    }
    std::vector<T> res(n, 0);
    matrix.multiply(v, res);
    ASSERT_EQ(res[0], -10);
    ASSERT_EQ(res[1], 8);
    ASSERT_EQ(res[2], 26);
}

TEST(Matrix, SquareComplexMultiplication){
    typedef std::complex<double> T;
    const size_t n=3;
    std::vector<T> v(n, 0);
    Matrix<T> matrix(n);
    for(size_t i=0; i<n; ++i) {
        v[i] = T{double(i+1), 2};
        for(size_t j=0; j<n; ++j) {
            matrix(i, j) = T{i*3.0-j*2.0+1, double(i+j)};
        }
    }
    std::vector<T> res(n, 0);
    matrix.multiply(v, res);
    ASSERT_EQ(res[0].real(), -16);
    ASSERT_EQ(res[0].imag(), 2);
    ASSERT_EQ(res[1].real(), -4);
    ASSERT_EQ(res[1].imag(), 26);
    ASSERT_EQ(res[2].real(), 8);
    ASSERT_EQ(res[2].imag(), 50);
}