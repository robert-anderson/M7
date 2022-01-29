//
// Created by rja on 16/06/2020.
//

#include <gtest/gtest.h>
#include "src/core/linalg/Dense.h"

TEST(Dense, Transpose) {
    dense::Matrix<int> p(3, 5);
    p.set_row(0, {1, 2, 3, 4, 5});
    p.set_row(1, {4, 2, 2, 7, 1});
    p.set_row(2, {8, 3, 4, 2, 7});
    auto p_chk = p;
    ASSERT_EQ(p_chk.nrow(), 3);
    ASSERT_EQ(p_chk.ncol(), 5);
    p.transpose();
    ASSERT_EQ(p(2, 1), 2);
    ASSERT_EQ(p(0, 2), 8);
    ASSERT_EQ(p.nrow(), 5);
    ASSERT_EQ(p.ncol(), 3);
    p.transpose();
    ASSERT_EQ(p.nrow(), 3);
    ASSERT_EQ(p.ncol(), 5);
    // two transposes should get back to original matrix
    ASSERT_TRUE(p.nearly_equal(p_chk));
}


TEST(Dense, Dagger) {
    dense::Matrix<std::complex<float>> p(3, 5);

    p.set_row(0, {{6.6, 8.4}, {8.7, 4.6}, {9.8, 7.8}, {4.7, 6.4}, {3.5, 6.0}});
    p.set_row(1, {{2.8, 3.2}, {5.6, 3.5}, {5.4, 1.5}, {5.8, 5.9}, {4.1, 5.4}});
    p.set_row(2, {{0.7, 6.4}, {4.8, 7.8}, {4.9, 7.0}, {7.4, 7.7}, {1.0, 7.1}});
    auto p_chk = p;
    ASSERT_EQ(p_chk.nrow(), 3);
    ASSERT_EQ(p_chk.ncol(), 5);
    p.dagger();
    ASSERT_FLOAT_EQ(p(2, 1).real(), 5.4);
    ASSERT_FLOAT_EQ(p(2, 1).imag(), -1.5);
    ASSERT_FLOAT_EQ(p(0, 2).real(), 0.7);
    ASSERT_FLOAT_EQ(p(0, 2).imag(), -6.4);
    ASSERT_EQ(p.nrow(), 5);
    ASSERT_EQ(p.ncol(), 3);
    p.dagger();
    ASSERT_EQ(p.nrow(), 3);
    ASSERT_EQ(p.ncol(), 5);
    // two transposes should get back to original matrix
    ASSERT_TRUE(p.nearly_equal(p_chk));
}

TEST(Dense, RectSingleMatMatProduct) {
    const float tol = 1e-6;
    dense::Matrix<float> p(3, 5);

    p.set_row(0, {0.7, 3. , 2.8, 2.9, 3.7});
    p.set_row(1, {0.6, 0.8, 0. , 2.6, 4. });
    p.set_row(2, {1.6, 2.2, 0.2, 3.3, 2.3});

    dense::Matrix<float> q(5, 4);

    q.set_row(0, {3.7, 2.7, 1.4, 1.4});
    q.set_row(1, {1.4, 3.3, 1.1, 1.2});
    q.set_row(2, {1.2, 3.4, 1.7, 3.8});
    q.set_row(3, {1. , 2.9, 2.6, 2.8});
    q.set_row(4, {3.4, 3.7, 1.1, 4. });

    dense::Matrix<float> r_chk(3, 4);
    r_chk.set_row(0, {25.63, 43.41, 20.65, 38.14});
    r_chk.set_row(1, {19.54, 26.6 , 12.88, 25.08});
    r_chk.set_row(2, {20.36, 30.34, 16.11, 24.08});

    dense::Matrix<float> r(3, 4);
    dense::multiply(p, q, r);
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));

    /*
     * check transpose flags in GEMM interface are working with P physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with P and Q physically transposed
     */
    r.zero();
    q.transpose();
    dense::multiply(p, q, r, 'T', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with Q physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'N', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
}

TEST(Dense, RectDoubleMatMatProduct) {
    const float tol = 1e-12;
    dense::Matrix<double> p(3, 5);

    p.set_row(0, {0.7, 3. , 2.8, 2.9, 3.7});
    p.set_row(1, {0.6, 0.8, 0. , 2.6, 4. });
    p.set_row(2, {1.6, 2.2, 0.2, 3.3, 2.3});

    dense::Matrix<double> q(5, 4);

    q.set_row(0, {3.7, 2.7, 1.4, 1.4});
    q.set_row(1, {1.4, 3.3, 1.1, 1.2});
    q.set_row(2, {1.2, 3.4, 1.7, 3.8});
    q.set_row(3, {1. , 2.9, 2.6, 2.8});
    q.set_row(4, {3.4, 3.7, 1.1, 4. });

    dense::Matrix<double> r_chk(3, 4);
    r_chk.set_row(0, {25.63, 43.41, 20.65, 38.14});
    r_chk.set_row(1, {19.54, 26.6 , 12.88, 25.08});
    r_chk.set_row(2, {20.36, 30.34, 16.11, 24.08});

    dense::Matrix<double> r(3, 4);
    dense::multiply(p, q, r);
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));

    /*
     * check transpose flags in GEMM interface are working with P physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with P and Q physically transposed
     */
    r.zero();
    q.transpose();
    dense::multiply(p, q, r, 'T', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with Q physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'N', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
}

TEST(Dense, RectComplexSingleMatMatProduct) {
    typedef float T;
    const T tol = 1e-5;
    dense::Matrix<std::complex<T>> p(3, 5);

    p.set_row(0, {{0.7, 1.6}, {3.0, 0.1}, {2.8, 3.0}, {2.9, 1.9}, {3.7, 3.6}});
    p.set_row(1, {{0.6, 4.0}, {0.8, 2.9}, {0.0, 3.5}, {2.6, 1.0}, {4.0, 1.8}});
    p.set_row(2, {{1.6, 2.0}, {2.2, 2.7}, {0.2, 0.2}, {3.3, 2.5}, {2.3, 2.0}});

    dense::Matrix<std::complex<T>> q(5, 4);

    q.set_row(0, {{0.5, 2.0}, {1.3, 0.7}, {1.5, 3.8}, {1.3, 0.9}});
    q.set_row(1, {{1.2, 1.5}, {1.5, 2.5}, {1.2, 2.3}, {0.1, 1.4}});
    q.set_row(2, {{3.8, 3.0}, {1.5, 3.3}, {0.9, 2.9}, {0.1, 2.9}});
    q.set_row(3, {{2.1, 1.4}, {2.9, 1.2}, {0.8, 1.1}, {1.2, 1.7}});
    q.set_row(4, {{2.0, 0.5}, {2.0, 2.2}, {1.5, 1.6}, {3.5, 0.4}});

    dense::Matrix<std::complex<T>> r_chk(3, 4);

    r_chk.set_row(0, {{11.27, 43.72}, {3.95, 48.29}, {-7.82, 38.93}, {2.97, 36.63}});
    r_chk.set_row(1, {{-10.43, 32.52}, {-9.24, 35.64}, {-26.06, 29.51}, {-2.25, 21.02}});
    r_chk.set_row(2, {{2.58, 27.12}, {3.64, 34.5} , {-9.03, 30.45}, {3.12, 24.52}});

    dense::Matrix<std::complex<T>> r(3, 4);
    dense::multiply(p, q, r);
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));

    /*
     * check transpose flags in GEMM interface are working with P physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with P and Q physically transposed
     */
    r.zero();
    q.transpose();
    dense::multiply(p, q, r, 'T', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with Q physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'N', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with P and Q physically conjugate-transposed
     */
    r.zero();
    p.dagger();
    q.conj();
    dense::multiply(p, q, r, 'C', 'C');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
}

TEST(Dense, RectComplexDoubleMatMatProduct) {
    typedef double T;
    const T tol = 1e-11;
    dense::Matrix<std::complex<T>> p(3, 5);

    p.set_row(0, {{0.7, 1.6}, {3.0, 0.1}, {2.8, 3.0}, {2.9, 1.9}, {3.7, 3.6}});
    p.set_row(1, {{0.6, 4.0}, {0.8, 2.9}, {0.0, 3.5}, {2.6, 1.0}, {4.0, 1.8}});
    p.set_row(2, {{1.6, 2.0}, {2.2, 2.7}, {0.2, 0.2}, {3.3, 2.5}, {2.3, 2.0}});

    dense::Matrix<std::complex<T>> q(5, 4);

    q.set_row(0, {{0.5, 2.0}, {1.3, 0.7}, {1.5, 3.8}, {1.3, 0.9}});
    q.set_row(1, {{1.2, 1.5}, {1.5, 2.5}, {1.2, 2.3}, {0.1, 1.4}});
    q.set_row(2, {{3.8, 3.0}, {1.5, 3.3}, {0.9, 2.9}, {0.1, 2.9}});
    q.set_row(3, {{2.1, 1.4}, {2.9, 1.2}, {0.8, 1.1}, {1.2, 1.7}});
    q.set_row(4, {{2.0, 0.5}, {2.0, 2.2}, {1.5, 1.6}, {3.5, 0.4}});

    dense::Matrix<std::complex<T>> r_chk(3, 4);

    r_chk.set_row(0, {{11.27, 43.72}, {3.95, 48.29}, {-7.82, 38.93}, {2.97, 36.63}});
    r_chk.set_row(1, {{-10.43, 32.52}, {-9.24, 35.64}, {-26.06, 29.51}, {-2.25, 21.02}});
    r_chk.set_row(2, {{2.58, 27.12}, {3.64, 34.5} , {-9.03, 30.45}, {3.12, 24.52}});

    dense::Matrix<std::complex<T>> r(3, 4);
    dense::multiply(p, q, r);
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));

    /*
     * check transpose flags in GEMM interface are working with P physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with P and Q physically transposed
     */
    r.zero();
    q.transpose();
    dense::multiply(p, q, r, 'T', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with Q physically transposed
     */
    r.zero();
    p.transpose();
    dense::multiply(p, q, r, 'N', 'T');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
    /*
     * check transpose flags in GEMM interface are working with P and Q physically conjugate-transposed
     */
    r.zero();
    p.dagger();
    q.conj();
    dense::multiply(p, q, r, 'C', 'C');
    ASSERT_TRUE(r.nearly_equal(r_chk, tol));
}

TEST(Dense, RectSingleMatVecProduct) {
    typedef float T;
    const T tol = 1e-5;
    dense::Matrix<T> p(3, 5);

    p.set_row(0, {0.7, 3. , 2.8, 2.9, 3.7});
    p.set_row(1, {0.6, 0.8, 0. , 2.6, 4. });
    p.set_row(2, {1.6, 2.2, 0.2, 3.3, 2.3});

    std::vector<T> q = {1.4, 3.3, 1.1, 1.2, 4.0};
    std::vector<T> r_chk = {32.24, 22.6 , 22.88};
    auto r = dense::multiply(p, q);
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));

    /*
     * now premultiply the vector by the transpose of the matrix
     */
    q = {1.4, 3.3, 1.1};
    r_chk = {4.72,  9.26,  4.14, 16.27, 20.91};
    r = dense::multiply(p, q, 'T');
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));
}

TEST(Dense, RectDoubleMatVecProduct) {
    typedef double T;
    const T tol = 1e-11;
    dense::Matrix<T> p(3, 5);

    p.set_row(0, {0.7, 3. , 2.8, 2.9, 3.7});
    p.set_row(1, {0.6, 0.8, 0. , 2.6, 4. });
    p.set_row(2, {1.6, 2.2, 0.2, 3.3, 2.3});

    std::vector<T> q = {1.4, 3.3, 1.1, 1.2, 4.0};
    std::vector<T> r_chk = {32.24, 22.6 , 22.88};
    auto r = dense::multiply(p, q);
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));

    /*
     * now premultiply a vector by the transpose of the matrix
     */
    q = {1.4, 3.3, 1.1};
    r_chk = {4.72,  9.26,  4.14, 16.27, 20.91};
    r = dense::multiply(p, q, 'T');
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));
}


TEST(Dense, RectComplexSingleMatVecProduct) {
    typedef float T;
    const T tol = 1e-5;
    dense::Matrix<std::complex<T>> p(3, 5);

    p.set_row(0, {{0.7, 1.6}, {3.0, 0.1}, {2.8, 3.0}, {2.9, 1.9}, {3.7, 3.6}});
    p.set_row(1, {{0.6, 4.0}, {0.8, 2.9}, {0.0, 3.5}, {2.6, 1.0}, {4.0, 1.8}});
    p.set_row(2, {{1.6, 2.0}, {2.2, 2.7}, {0.2, 0.2}, {3.3, 2.5}, {2.3, 2.0}});

    std::vector<std::complex<T>> q = {{2.6, 2.5}, {3.9, 1.3}, {2.5, 4.0}, {2.9, 2.2}, {2.1, 1.4}};
    std::vector<std::complex<T>> r_chk = {{11.35, 53.53}, {-11.87, 51.0}, {10.03, 45.82}};
    auto r = dense::multiply(p, q);
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));

    /*
     * now premultiply a vector by the transpose of the matrix
     */
    q = {{2.6, 2.5}, {3.9, 1.3}, {2.5, 4.0}};
    r_chk = {{-9.04, 33.69}, {1.6, 35.66}, {-5.07, 29.75}, {9.88, 38.92}, {11.63, 45.03}};
    r = dense::multiply(p, q, 'T');
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));

    /*
     * finally, premultiply the same vector by the conjugate-transpose of the matrix
     */
    r_chk = {{25.36, -15.83}, {31.24, -0.98}, {20.63, -14.15}, {41.98, 8.74}, {50.31, 2.27}};
    r = dense::multiply(p, q, 'C');
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));
}

TEST(Dense, RectComplexDoubleMatVecProduct) {
    typedef double T;
    const T tol = 1e-11;
    dense::Matrix<std::complex<T>> p(3, 5);

    p.set_row(0, {{0.7, 1.6}, {3.0, 0.1}, {2.8, 3.0}, {2.9, 1.9}, {3.7, 3.6}});
    p.set_row(1, {{0.6, 4.0}, {0.8, 2.9}, {0.0, 3.5}, {2.6, 1.0}, {4.0, 1.8}});
    p.set_row(2, {{1.6, 2.0}, {2.2, 2.7}, {0.2, 0.2}, {3.3, 2.5}, {2.3, 2.0}});

    std::vector<std::complex<T>> q = {{2.6, 2.5}, {3.9, 1.3}, {2.5, 4.0}, {2.9, 2.2}, {2.1, 1.4}};
    std::vector<std::complex<T>> r_chk = {{11.35, 53.53}, {-11.87, 51.0}, {10.03, 45.82}};
    auto r = dense::multiply(p, q);
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));

    /*
     * now premultiply a vector by the transpose of the matrix
     */
    q = {{2.6, 2.5}, {3.9, 1.3}, {2.5, 4.0}};
    r_chk = {{-9.04, 33.69}, {1.6, 35.66}, {-5.07, 29.75}, {9.88, 38.92}, {11.63, 45.03}};
    r = dense::multiply(p, q, 'T');
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));

    /*
     * finally, premultiply the same vector by the conjugate-transpose of the matrix
     */
    r_chk = {{25.36, -15.83}, {31.24, -0.98}, {20.63, -14.15}, {41.98, 8.74}, {50.31, 2.27}};
    r = dense::multiply(p, q, 'C');
    ASSERT_TRUE(dense::nearly_equal(r_chk, r, tol));
}

TEST(Dense, SingleInnerProduct) {
    typedef float T;
    const T tol = 1e-5;
    std::vector<T> p = {0.7, 3. , 2.8, 2.9, 3.7};
    std::vector<T> q = {0.6, 0.8, 0. , 2.6, 4. };
    ASSERT_TRUE(consts::nearly_equal(dense::inner_product(p, q), T(25.16), tol));
}

TEST(Dense, DoubleInnerProduct) {
    typedef double T;
    const T tol = 1e-11;
    std::vector<T> p = {0.7, 3. , 2.8, 2.9, 3.7};
    std::vector<T> q = {0.6, 0.8, 0. , 2.6, 4. };
    ASSERT_TRUE(consts::nearly_equal(dense::inner_product(p, q), T(25.16), tol));
}

TEST(Dense, ComplexSingleInnerProduct) {
    typedef float T;
    const T tol = 1e-5;
    std::vector<std::complex<T>> p = {{0.7, 1.6}, {3.0, 0.1}, {2.8, 3.0}, {2.9, 1.9}, {3.7, 3.6}};
    std::vector<std::complex<T>> q = {{0.6, 4.0}, {0.8, 2.9}, {0.0, 3.5}, {2.6, 1.0}, {4.0, 1.8}};
    auto r = dense::inner_product(p, q);
    ASSERT_TRUE(consts::nearly_equal(r.real(), T(-0.41), tol));
    ASSERT_TRUE(consts::nearly_equal(r.imag(), T(51.24), tol));
    r = dense::inner_product(p, q, true);
    ASSERT_TRUE(consts::nearly_equal(r.real(), T(50.73), tol));
    ASSERT_TRUE(consts::nearly_equal(r.imag(), T(10.48), tol));
}

TEST(Dense, ComplexDoubleInnerProduct) {
    typedef double T;
    const T tol = 1e-11;
    std::vector<std::complex<T>> p = {{0.7, 1.6}, {3.0, 0.1}, {2.8, 3.0}, {2.9, 1.9}, {3.7, 3.6}};
    std::vector<std::complex<T>> q = {{0.6, 4.0}, {0.8, 2.9}, {0.0, 3.5}, {2.6, 1.0}, {4.0, 1.8}};
    auto r = dense::inner_product(p, q);
    ASSERT_TRUE(consts::nearly_equal(r.real(), T(-0.41), tol));
    ASSERT_TRUE(consts::nearly_equal(r.imag(), T(51.24), tol));
    r = dense::inner_product(p, q, true);
    ASSERT_TRUE(consts::nearly_equal(r.real(), T(50.73), tol));
    ASSERT_TRUE(consts::nearly_equal(r.imag(), T(10.48), tol));
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

#if 0

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

#endif
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