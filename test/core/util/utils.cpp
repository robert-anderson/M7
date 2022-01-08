//
// Created by rja on 07/08/2021.
//

#include "src/core/util/utils.h"
#include "gtest/gtest.h"

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

TEST(utils, SplitLine){
    const std::string line = "0.5000000000 1 1 2 2";
    auto tokens = string_utils::split(line, ' ');
    std::cout << tokens << std::endl;
}

TEST(utils, JoinAndSplit){
    const std::string line = " this is   an   example   string   ";
    auto tokens = string_utils::split(line, ' ');
    ASSERT_EQ(tokens.size(), 5);
    auto joinder = string_utils::join(tokens, " ");
    // splitting will eliminate consecutive occurrences of the delimiter
    ASSERT_EQ("this is an example string", joinder);
}

TEST(utils, Tokenize){
    const std::string line = " this is   an,   example   string   ";
    auto tokens = string_utils::split(line, " ,");
    ASSERT_EQ(tokens.size(), 5);
}

TEST(Utils, CompileTimePow){
    ASSERT_EQ(utils::pow<3>(5), 5*5*5);
    ASSERT_EQ(utils::pow<10>(2), 1024);
    ASSERT_EQ(utils::pow<0>(10), 1ul);
    ASSERT_EQ(utils::pow<1>(10), 10ul);
    ASSERT_EQ(utils::pow<1>(0), 0ul);
}
TEST(Utils, CompileTimeNtup){
    ASSERT_EQ(utils::ntup<4>(15), integer_utils::combinatorial(15, 4));
    ASSERT_EQ(utils::ntup<1>(15), 15);
    ASSERT_EQ(utils::ntup<1>(1), 1);
    ASSERT_EQ(utils::ntup<0>(1), 1);
    ASSERT_EQ(utils::ntup<0>(15), 1);
}

TEST(Utils, Exsigs) {
    using namespace exsig_utils;
    /*
     * assert that all excitation signatures are encoded and decoded correctly
     */
    size_t exsig;
    for (size_t ncref = 0ul; ncref<=defs::exsig_nop_mask_frm; ++ncref){
        for (size_t nannf = 0ul; nannf<=defs::exsig_nop_mask_frm; ++nannf){
            for (size_t ncreb = 0ul; ncreb<=defs::exsig_nop_mask_bos; ++ncreb){
                for (size_t nannb = 0ul; nannb<=defs::exsig_nop_mask_bos; ++nannb){
                    exsig = encode(ncref, nannf, ncreb, nannb);
                    ASSERT_EQ(decode_nfrm_cre(exsig), ncref);
                    ASSERT_EQ(decode_nfrm_ann(exsig), nannf);
                    ASSERT_EQ(decode_nbos_cre(exsig), ncreb);
                    ASSERT_EQ(decode_nbos_ann(exsig), nannb);
                }
                exsig = encode(ncref, nannf, ncreb, defs::exsig_nop_mask_bos + 1);
                ASSERT_EQ(exsig, ~0ul);
            }
            exsig = encode(ncref, nannf, defs::exsig_nop_mask_bos + 1, 0ul);
            ASSERT_EQ(exsig, ~0ul);
        }
        exsig = encode(ncref, defs::exsig_nop_mask_frm + 1, 0ul, 0ul);
        ASSERT_EQ(exsig, ~0ul);
    }
    exsig = encode(defs::exsig_nop_mask_frm + 1, 0ul, 0ul, 0ul);
    ASSERT_EQ(exsig, ~0ul);
}


#if 0
TEST(Utils, SetAllExsigsFromRanksig) {
    size_t ranksig;
    std::array<bool, defs::nexsig> exsigs{};

    ranksig = conn_utils::encode_exsig(4, 4, 1, 1);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(4, 4, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(4, 4, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 0, 0)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 10);

    ranksig = conn_utils::encode_exsig(4, 4, 0, 0);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(4, 4, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 0, 0)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 5);

    ranksig = conn_utils::encode_exsig(3, 3, 1, 0);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 1, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 1, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 1, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 1, 0)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 4);


    ranksig = conn_utils::encode_exsig(2, 2, 0, 1);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 0, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 0, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 0, 1)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 3);
}
#endif