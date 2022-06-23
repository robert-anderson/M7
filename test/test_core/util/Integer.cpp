
//
// Created by Robert J. Anderson on 07/08/2021.
//

#include "test_core/defs.h"
#include "M7_lib/util/Integer.h"
#include "M7_lib/util/MeanStd.h"

using namespace integer;

TEST(UtilInteger, Factorial) {
    ASSERT_EQ(factorial(0), 1ul);
    ASSERT_EQ(factorial(1), 1ul);
    ASSERT_EQ(factorial(2), 2ul);
    ASSERT_EQ(factorial(3), 6ul);
    ASSERT_EQ(factorial(20), 2432902008176640000ul);
    // overflows at 21!
    ASSERT_EQ(factorial(21), 14197454024290336768ul);
}

TEST(UtilInteger, Combinatorial) {
    for (size_t n = 0ul; n < 20; ++n) {
        for (size_t r = 0ul; r <= n; ++r) {
            ASSERT_EQ(combinatorial(n, r),
                      factorial(n) / (factorial(n - r) * factorial(r)));
        }
    }
    ASSERT_EQ(combinatorial(100, 5), 75287520ul);
    ASSERT_EQ(combinatorial(100, 10), 17310309456440ul);
}

TEST(UtilInteger, PairMaps) {
    const size_t N = 20;
    size_t n;
    size_t tmp_i, tmp_j;

    size_t ij = 0ul;
    for (size_t i = 0ul; i < N; ++i) {
        for (size_t j = 0ul; j <= i; ++j) {
            n = trigmap(i, j);
            inv_trigmap(tmp_i, tmp_j, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }

    ij = 0ul;
    for (size_t i = 0ul; i < N; ++i) {
        for (size_t j = 0ul; j < i; ++j) {
            n = strigmap(i, j);
            inv_strigmap(tmp_i, tmp_j, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }

    ij = 0ul;
    for (size_t i = 0ul; i < N; ++i) {
        for (size_t j = 0ul; j < N; ++j) {
            n = rectmap(i, j, N);
            inv_rectmap(tmp_i, tmp_j, N, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }
}

TEST(UtilInteger, CompileTimeNtup) {
    ASSERT_EQ(ntup<4>(15), integer::combinatorial(15, 4));
    ASSERT_EQ(ntup<1>(15), 15);
    ASSERT_EQ(ntup<1>(1), 1);
    ASSERT_EQ(ntup<0>(1), 1);
    ASSERT_EQ(ntup<0>(15), 1);
}