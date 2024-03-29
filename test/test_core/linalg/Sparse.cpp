//
// Created by Robert J. Anderson on 07/01/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/linalg/sparse/Inverse.h"

TEST(Sparse, RectMatrixVectorProduct){
    typedef int T;
    sparse::dynamic::Matrix<T> mat;
    /*
     *  a[0,0], a[0,1], a[1,2], a[2,2], a[3,0], a[4,1], a[4,2], a[5, 0], a[5,2], a[6, 1] =
     *      1, -1, 2, -2, 3, -3, 4, -4, 5, -5
     *  v = [-9, 8, 2]
     *
     *  a.v = [-17, 4, -4, -27, -16,  46, -40]
     */
    const uint_t nrow = 7;
    uintv_t irows = {0, 0, 1, 2, 3, 4, 4, 5, 5, 6};
    uintv_t icols = {0, 1, 2, 2, 0, 1, 2, 0, 2, 1};
    v_t<T> vs = {1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
    mat.resize(nrow);
    for (uint_t i=0ul; i<irows.size(); ++i) mat.add(irows[i], {icols[i], vs[i]});

    v_t<T> in = {-9, 8, 2};
    ASSERT_EQ(in.size(), mat.max_col_ind()+1);
    v_t<T> out;
    v_t<T> out_chk = {-17, 4, -4, -27, -16,  46, -40};
    mat.multiply(in, out);
    ASSERT_EQ(out, out_chk);
}

TEST(Sparse, DynamicToFixed){
    typedef float T;
    v_t<std::pair<uint_t, sparse::MatrixElement<T>>> contents {
            {11, {14, 0.12}},
            {11, {12, 1.23}},
            {11, {16, 2.34}},
            {11, {13, 3.45}},
            {11, {2, 4.56}},
            {3, {2, 5.67}},
            {3, {1, 6.78}},
            {3, {11, 7.89}}
    };

    v_t<std::pair<uint_t, sparse::MatrixElement<T>>> contents_sorted {
            {11, {2, 4.56}},
            {11, {12, 1.23}},
            {11, {13, 3.45}},
            {11, {14, 0.12}},
            {11, {16, 2.34}},
            {3, {1, 6.78}},
            {3, {2, 5.67}},
            {3, {11, 7.89}}
    };

    sparse::dynamic::Matrix<T> dynamic;
    for (auto& elem: contents) dynamic.insert(elem.first, elem.second);
    sparse::dynamic::Matrix<T> dynamic_sorted;
    for (auto& elem: contents_sorted) dynamic_sorted.insert(elem.first, elem.second);

    sparse::fixed::Matrix<T> fixed(dynamic_sorted, false);
    sparse::fixed::Matrix<T> fixed_sorted(dynamic, true);

    ASSERT_EQ(fixed, fixed_sorted);
}