//
// Created by Robert John Anderson on 2020-01-18.
//

#include <gtest/gtest.h>
#include "src/core/util/utils.h"

using namespace integer_utils;

TEST(utils, factorial){
    ASSERT_EQ(factorial(0), 1ul);
    ASSERT_EQ(factorial(1), 1ul);
    ASSERT_EQ(factorial(2), 2ul);
    ASSERT_EQ(factorial(3), 6ul);
    ASSERT_EQ(factorial(20), 2432902008176640000ul);
    // overflows at 21!
    ASSERT_EQ(factorial(21), 14197454024290336768ul);
}

TEST(utils, combinatorial){
    for (size_t n = 0ul; n<20; ++n) {
        for (size_t r = 0ul; r <= n; ++r) {
            ASSERT_EQ(combinatorial(n, r),
                      factorial(n) / (factorial(n - r) * factorial(r)));
        }
    }
    ASSERT_EQ(combinatorial(100, 5), 75287520ul);
    ASSERT_EQ(combinatorial(100, 10), 17310309456440ul);
}

TEST(utils, PairMaps){
    const size_t N = 20;
    size_t n;
    size_t tmp_i, tmp_j;

    size_t ij = 0ul;
    for (size_t i=0ul; i<N; ++i){
        for(size_t j=0ul; j<=i; ++j){
            n = trigmap(i, j);
            inv_trigmap(tmp_i, tmp_j, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }

    ij = 0ul;
    for (size_t i=0ul; i<N; ++i){
        for(size_t j=0ul; j<i; ++j){
            n = strigmap(i, j);
            inv_strigmap(tmp_i, tmp_j, n);
            ASSERT_EQ(n, ij);
            ASSERT_EQ(tmp_i, i);
            ASSERT_EQ(tmp_j, j);
            ++ij;
        }
    }

    ij = 0ul;
    for (size_t i=0ul; i<N; ++i) {
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

TEST(utils, MeanAndStd){
    std::vector<double> v = {1, 2, 3.8, 4, -0.35, 0.6};
    auto mean_std = stat_utils::mean_std<double>(v.cbegin(), v.cend());
    ASSERT_FLOAT_EQ(mean_std.first,  1.8416666666666668);
    ASSERT_FLOAT_EQ(mean_std.second, 1.6110081384717527);
}

TEST(utils, JoinAndSplit){
    const std::string line = " this is   an   example   string   ";
    auto tokens = string_utils::split(line, ' ');
    ASSERT_EQ(tokens.size(), 5);
    auto joinder = string_utils::join(tokens, " ", false);
    std::cout << joinder << std::endl;
    // splitting will eliminate consecutive occurrences of the delimiter
    ASSERT_EQ("this is an example string", joinder);
}

TEST(utils, Tokenize){
    const std::string line = " this is   an,   example   string   ";
    auto tokens = string_utils::split(line, " ,");
    ASSERT_EQ(tokens.size(), 5);
}
