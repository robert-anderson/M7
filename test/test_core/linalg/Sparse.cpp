//
// Created by Robert J. Anderson on 07/01/2022.
//

#include "M7_lib/linalg/Sparse.h"
#include "gtest/gtest.h"

TEST(Sparse, RectMatrixVectorProduct){
    typedef int T;
    sparse::Matrix<T> mat;
    /*
     *  a[0,0], a[0,1], a[1,2], a[2,2], a[3,0], a[4,1], a[4,2], a[5, 0], a[5,2], a[6, 1] =
     *      1, -1, 2, -2, 3, -3, 4, -4, 5, -5
     *  v = [-9, 8, 2]
     *
     *  a.v = [-17, 4, -4, -27, -16,  46, -40]
     */
    const size_t nrow = 7;
    defs::inds irows = {0, 0, 1, 2, 3, 4, 4, 5, 5, 6};
    defs::inds icols = {0, 1, 2, 2, 0, 1, 2, 0, 2, 1};
    std::vector<T> vs = {1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
    mat.resize(nrow);
    for (size_t i=0ul; i<irows.size(); ++i) mat.add(irows[i], icols[i], vs[i]);

    std::vector<T> in = {-9, 8, 2};
    ASSERT_EQ(in.size(), mat.max_column_index()+1);
    std::vector<T> out;
    std::vector<T> out_chk = {-17, 4, -4, -27, -16,  46, -40};
    mat.multiply(in, out);
    ASSERT_EQ(out, out_chk);
}