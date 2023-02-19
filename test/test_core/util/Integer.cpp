
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
    for (uint_t n = 0ul; n < 20; ++n) {
        for (uint_t r = 0ul; r <= n; ++r) {
            ASSERT_EQ(combinatorial(n, r),
                      factorial(n) / (factorial(n - r) * factorial(r)));
        }
    }
    ASSERT_EQ(combinatorial(100, 5), 75287520ul);
    ASSERT_EQ(combinatorial(100, 10), 17310309456440ul);
}

TEST(UtilInteger, PairMaps) {
    const uint_t N = 20;
    uint_t n;
    uint_t tmp_i, tmp_j;

    uint_t ij = 0ul;
    for (uint_t i = 0ul; i < N; ++i) {
        for (uint_t j = 0ul; j <= i; ++j) {
            n = trigmap(i, j);
            inv_trigmap(tmp_i, tmp_j, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }

    ij = 0ul;
    for (uint_t i = 0ul; i < N; ++i) {
        for (uint_t j = 0ul; j < i; ++j) {
            n = strigmap(i, j);
            inv_strigmap(tmp_i, tmp_j, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }

    ij = 0ul;
    for (uint_t i = 0ul; i < N; ++i) {
        for (uint_t j = 0ul; j < N; ++j) {
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

TEST(UtilInteger, LeastCommonMultipleLe) {
    ASSERT_EQ(integer::lcm_le(1), 1);
    ASSERT_EQ(integer::lcm_le(2), 2);
    ASSERT_EQ(integer::lcm_le(3), 3*2);
    ASSERT_EQ(integer::lcm_le(4), 3*2*2);
    ASSERT_EQ(integer::lcm_le(5), 5*3*2*2);
    ASSERT_EQ(integer::lcm_le(6), 5*3*2*2);
    ASSERT_EQ(integer::lcm_le(7), 7*5*3*2*2);
    ASSERT_EQ(integer::lcm_le(8), 7*5*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(9), 7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(10), 7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(11), 11*7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(12), 11*7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(13), 13*11*7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(14), 13*11*7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(15), 13*11*7*5*3*3*2*2*2);
    ASSERT_EQ(integer::lcm_le(16), 13*11*7*5*3*3*2*2*2*2);
}

TEST(UtilInteger, Sqrt) {
    const uint_t nroot = 100ul;
    uint_t inonsq = 0ul;
    for (uint_t iroot = 0ul; iroot<nroot; ++iroot){
        const auto isq = iroot * iroot;
        for (; inonsq != isq; ++inonsq){
            ASSERT_EQ(integer::sqrt(inonsq), ~0ul);
        }
        ++inonsq;
        ASSERT_EQ(integer::sqrt(isq), iroot);
    }
}

TEST(UtilInteger, Partition) {
    {
        auto parts = partitions(1);
        v_t<uintv_t> chk = {
            {1}
        };
        ASSERT_EQ(parts, chk);
    }
    {
        auto parts = partitions(2);
        v_t<uintv_t> chk = {
            {2},
            {1, 1}
        };
        ASSERT_EQ(parts, chk);
    }
    {
        auto parts = partitions(3);
        v_t<uintv_t> chk = {
            {3},
            {2, 1},
            {1, 1, 1}
        };
        ASSERT_EQ(parts, chk);
    }
    {
        auto parts = partitions(6);
        v_t<uintv_t> chk = {
            {6},
            {5, 1},
            {4, 2},
            {4, 1, 1},
            {3, 3},
            {3, 2, 1},
            {3, 1, 1, 1},
            {2, 2, 2},
            {2, 2, 1, 1},
            {2, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1}
        };
        ASSERT_EQ(parts, chk);
    }
    {
        auto parts = partitions(6, 3);
        v_t<uintv_t> chk = {
            {3, 3},
            {3, 2, 1},
            {3, 1, 1, 1},
            {2, 2, 2},
            {2, 2, 1, 1},
            {2, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1}
        };
        ASSERT_EQ(parts, chk);
    }
}