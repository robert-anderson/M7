//
// Created by rja on 16/06/2020.
//

#include <gtest/gtest.h>
#include "src/core/linalg/Dense.h"

TEST(Dense, SquareRealMultiplication) {
    typedef double T;
    const size_t n = 3;
    std::vector<T> v(n, 0);
    dense::SquareMatrix<T> matrix(n);
    for (size_t i = 0; i < n; ++i) {
        v[i] = i + 1;
        for (size_t j = 0; j < n; ++j) {
            matrix(i, j) = i * 3.0 - j * 2.0 + 1;
        }
    }
    std::vector<T> res(n, 0);
    matrix.multiply(v, res);
    ASSERT_EQ(res[0], -10);
    ASSERT_EQ(res[1], 8);
    ASSERT_EQ(res[2], 26);
}

TEST(Dense, SquareComplexMultiplication) {
    typedef std::complex<double> T;
    const size_t n = 3;
    std::vector<T> v(n, 0);
    dense::SquareMatrix<T> matrix(n);
    for (size_t i = 0; i < n; ++i) {
        v[i] = T{double(i + 1), 2};
        for (size_t j = 0; j < n; ++j) {
            matrix(i, j) = T{i * 3.0 - j * 2.0 + 1, double(i + j)};
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

TEST(Dense, RealSymEig) {
    typedef double T;
    const size_t n = 4;
    dense::SquareMatrix<T> mat(n);
    /*
     * [[ 0,  2, -1,  2],
       [ 2,  4, -1, -1],
       [-1, -1,  2, -5],
       [ 2, -1, -5,  4]]
     */
    mat.set_row(0, {0, 2, -1, 2});
    mat.set_row(1, {2, 4, -1, -1});
    mat.set_row(2, {-1, -1, 2, -5});
    mat.set_row(3, {2, -1, -5, 4});

    auto evals = mat.diag_inplace<double>();

    std::vector<T> evals_chk = {-2.84695223, -0.7346574, 4.90333541, 8.67827421};

    for (size_t i = 0ul; i < n; ++i) ASSERT_FLOAT_EQ(evals[i], evals_chk[i]);
    /*
     * [[ 0.42674891, -0.30080954, -0.59964891, -0.6064818 ],
       [ 0.81540297, -0.2009709 ,  0.51855174,  0.16072585],
       [-0.29154909, -0.92861347,  0.02794618,  0.22780509],
       [ 0.26077289,  0.08247022, -0.60888775,  0.74461525]]
     */
    dense::SquareMatrix<T> evecs(n);
    evecs.set_row(0, {0.42674891, -0.30080954, -0.59964891, -0.6064818});
    evecs.set_row(1, {0.81540297, -0.2009709, 0.51855174, 0.16072585});
    evecs.set_row(2, {-0.29154909, -0.92861347, 0.02794618, 0.22780509});
    evecs.set_row(3, {0.26077289, 0.08247022, -0.60888775, 0.74461525});

    ASSERT_TRUE(mat.nearly_equal(evecs, 1e-7));
}


TEST(Dense, ComplexHermEig) {
    typedef double T;
    const size_t n = 4;
    dense::SquareMatrix<std::complex<T>> mat(n);
    /*
     * [[ 4.+0.j, -1.+0.j, -5.+2.j, -2.-5.j],
       [-1.+0.j,  0.+0.j, -3.-1.j, -3.+2.j],
       [-5.-2.j, -3.+1.j, -4.+0.j, -1.+2.j],
       [-2.+5.j, -3.-2.j, -1.-2.j,  2.+0.j]]
     */
    mat.set_row(0, {{4,  0},
                    {-1, 0},
                    {-5, 2},
                    {-2, -5}});
    mat.set_row(1, {{-1, 0},
                    {0,  0},
                    {-3, -1},
                    {-3, 2}});
    mat.set_row(2, {{-5, -2},
                    {-3, 1},
                    {-4, 0},
                    {-1, 2}});
    mat.set_row(3, {{-2, 5},
                    {-3, -2},
                    {-1, -2},
                    {2,  0}});

    auto evals = mat.diag_inplace<double>();

    std::vector<T> evals_chk = {-8.50478306, -3.4626037, 3.23302049, 10.73436627};

    for (size_t i = 0ul; i < n; ++i) ASSERT_FLOAT_EQ(evals[i], evals_chk[i]);

    /*
     * real:
     * [[ 3.58482045e-01,  5.10924587e-01,  2.72230923e-01,  -7.32350336e-01],
       [ 4.23494924e-01, -3.15364004e-01,  1.79580298e-01,    5.40389092e-02],
       [ 7.86848083e-01,  2.34349354e-02, -2.92244286e-01,    2.92874123e-01],
       [ 2.46562766e-01, -5.26100478e-04,  1.48523673e-01,    1.75533836e-01]]
     * imag
     * [[-0.        ,  0.        ,  0.        ,  0.        ],
       [ 0.07095225, -0.28298682,  0.76554493,  0.12187498],
       [ 0.02115204,  0.33401981, -0.27849933,  0.13985851],
       [ 0.08208495, -0.66880593, -0.3500776 , -0.55654424]]
     */
    std::vector<T> evecs_real = {
            3.58482045e-01, 5.10924587e-01, 2.72230923e-01, -7.32350336e-01,
            4.23494924e-01, -3.15364004e-01, 1.79580298e-01, 5.40389092e-02,
            7.86848083e-01, 2.34349354e-02, -2.92244286e-01, 2.92874123e-01,
            2.46562766e-01, -5.26100478e-04, 1.48523673e-01, 1.75533836e-01};
    std::vector<T> evecs_imag = {
            -0., 0., 0., 0.,
            0.07095225, -0.28298682, 0.76554493, 0.12187498,
            0.02115204, 0.33401981, -0.27849933, 0.13985851,
            0.08208495, -0.66880593, -0.3500776, -0.55654424};

    dense::SquareMatrix<std::complex<T>> evecs(n);
    for (size_t i=0ul; i<evecs_real.size(); ++i) evecs[i] = {evecs_real[i], evecs_imag[i]};

    ASSERT_TRUE(mat.nearly_equal(evecs, 1e-7));
}

#if 0
TEST(Dense, RealNonSymEig) {
    typedef double T;
    const size_t n=4;
    dense::SquareMatrix<T> mat(n);
    /*
     * [[-5,  0, -4,  0],
       [-3, -1,  1,  0],
       [-3,  0, -5, -2],
       [ 4, -2, -1, -5]]
     */
    mat.set_row(0, {-5,  0, -4,  0});
    mat.set_row(1, {-3, -1,  1,  0});
    mat.set_row(2, {-3,  0, -5, -2});
    mat.set_row(3, {4, -2, -1, -5});

    auto evals = mat.diag_inplace<double>();

    std::vector<std::complex<T>> evals_chk = {{-7.09850305, 0.94428128}, {-7.09850305, -0.94428128},
                                              0.39858951, -2.2015834};

    for (size_t i=0ul; i<n; ++i) {
        ASSERT_FLOAT_EQ(evals[i].real(), evals_chk[i].real());
        ASSERT_FLOAT_EQ(evals[i].imag(), evals_chk[i].imag());
    }
}


TEST(Dense, RealNonSymEigRight) {
    typedef double T;
    const size_t n=4;
    dense::SquareMatrix<T> mat(n);
    /*
     * [[-5,  0, -4,  0],
       [-3, -1,  1,  0],
       [-3,  0, -5, -2],
       [ 4, -2, -1, -5]]
     */
    mat.set_row(0, {-5,  0, -4,  0});
    mat.set_row(1, {-3, -1,  1,  0});
    mat.set_row(2, {-3,  0, -5, -2});
    mat.set_row(3, {4, -2, -1, -5});

    auto evals = mat.diag_inplace<std::complex<double>>();

    std::vector<std::complex<T>> evals_chk = {{-7.09850305, 0.94428128}, {-7.09850305, -0.94428128},
                                              0.39858951, -2.2015834};

    for (size_t i=0ul; i<n; ++i) {
        ASSERT_FLOAT_EQ(evals[i].real(), evals_chk[i].real());
        ASSERT_FLOAT_EQ(evals[i].imag(), evals_chk[i].imag());
    }
}
#endif