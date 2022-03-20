//
// Created by rja on 07/08/2021.
//

#include "src/core/util/utils.h"
#include "gtest/gtest.h"

using namespace integer_utils;

TEST(utils, Factorial) {
    ASSERT_EQ(factorial(0), 1ul);
    ASSERT_EQ(factorial(1), 1ul);
    ASSERT_EQ(factorial(2), 2ul);
    ASSERT_EQ(factorial(3), 6ul);
    ASSERT_EQ(factorial(20), 2432902008176640000ul);
    // overflows at 21!
    ASSERT_EQ(factorial(21), 14197454024290336768ul);
}

TEST(utils, Combinatorial) {
    for (size_t n = 0ul; n < 20; ++n) {
        for (size_t r = 0ul; r <= n; ++r) {
            ASSERT_EQ(combinatorial(n, r),
                      factorial(n) / (factorial(n - r) * factorial(r)));
        }
    }
    ASSERT_EQ(combinatorial(100, 5), 75287520ul);
    ASSERT_EQ(combinatorial(100, 10), 17310309456440ul);
}

TEST(utils, PairMaps) {
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

TEST(utils, MeanAndStd) {
    std::vector<double> v = {1, 2, 3.8, 4, -0.35, 0.6};
    auto mean_std = stat_utils::mean_std<double>(v.cbegin(), v.cend());
    ASSERT_FLOAT_EQ(mean_std.first, 1.8416666666666668);
    ASSERT_FLOAT_EQ(mean_std.second, 1.6110081384717527);
}

TEST(utils, SplitLine) {
    const std::string line = "0.5000000000 1 1 2 2";
    auto tokens = string_utils::split(line, ' ');
}

TEST(utils, JoinAndSplit) {
    const std::string line = " this is   an   example   string   ";
    auto tokens = string_utils::split(line, ' ');
    ASSERT_EQ(tokens.size(), 5);
    auto joinder = string_utils::join(tokens, " ");
    // splitting will eliminate consecutive occurrences of the delimiter
    ASSERT_EQ("this is an example string", joinder);
}

TEST(utils, Tokenize) {
    const std::string line = " this is   an,   example   string   ";
    auto tokens = string_utils::split(line, " ,");
    ASSERT_EQ(tokens.size(), 5);
}

TEST(Utils, CompileTimePow) {
    ASSERT_EQ(utils::pow<3>(5), 5 * 5 * 5);
    ASSERT_EQ(utils::pow<10>(2), 1024);
    ASSERT_EQ(utils::pow<0>(10), 1ul);
    ASSERT_EQ(utils::pow<1>(10), 10ul);
    ASSERT_EQ(utils::pow<1>(0), 0ul);
}

TEST(Utils, CompileTimeNtup) {
    ASSERT_EQ(utils::ntup<4>(15), integer_utils::combinatorial(15, 4));
    ASSERT_EQ(utils::ntup<1>(15), 15);
    ASSERT_EQ(utils::ntup<1>(1), 1);
    ASSERT_EQ(utils::ntup<0>(1), 1);
    ASSERT_EQ(utils::ntup<0>(15), 1);
}

TEST(Utils, SetBits) {
    defs::inds setbits;
    std::string bitstr;

    bitstr = "10011110"; // 158
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint8_t u8 = 0;
    for (const auto& bit: setbits) bit_utils::set(u8, bit);
    ASSERT_EQ(u8, 158);

    bitstr = "11110111000100"; // 15812
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint16_t u16 = 0;
    for (const auto& bit: setbits) bit_utils::set(u16, bit);
    ASSERT_EQ(u16, 15812);

    bitstr = "10100010111000001110110011101"; // 341581213
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint32_t u32 = 0;
    for (const auto& bit: setbits) bit_utils::set(u32, bit);
    ASSERT_EQ(u32, 341581213);

    bitstr = "10000001111110111101111011110000000001100111110101001100011110"; // 2341581243132367646
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint64_t u64 = 0;
    for (const auto& bit: setbits) bit_utils::set(u64, bit);
    ASSERT_EQ(u64, 2341581243132367646);
}


TEST(Utils, ClrBits) {
    defs::inds clrbits;
    std::string bitstr;

    bitstr = "10011110"; // 158
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint8_t u8 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit_utils::set(u8, ibit);
    for (const auto& bit: clrbits) bit_utils::clr(u8, bit);
    ASSERT_EQ(u8, 158);

    bitstr = "11110111000100"; // 15812
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint16_t u16 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit_utils::set(u16, ibit);
    for (const auto& bit: clrbits) bit_utils::clr(u16, bit);
    ASSERT_EQ(u16, 15812);

    bitstr = "10100010111000001110110011101"; // 341581213
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint32_t u32 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit_utils::set(u32, ibit);
    for (const auto& bit: clrbits) bit_utils::clr(u32, bit);
    ASSERT_EQ(u32, 341581213);

    bitstr = "10000001111110111101111011110000000001100111110101001100011110"; // 2341581243132367646
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint64_t u64 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit_utils::set(u64, ibit);
    for (const auto& bit: clrbits) bit_utils::clr(u64, bit);
    ASSERT_EQ(u64, 2341581243132367646);
}

TEST(Utils, BitRangeMasks) {
    uint8_t u8 = 0;
    uint16_t u16 = 0;
    uint32_t u32 = 0;
    uint64_t u64 = 0;
    /*
     * whole integer edge cases: negations should be zero
     */
    bit_utils::make_range_mask(u8, 0, 8);
    ASSERT_FALSE(uint8_t(~u8));
    bit_utils::make_range_mask(u16, 0, 16);
    ASSERT_FALSE(uint16_t(~u16));
    bit_utils::make_range_mask(u32, 0, 32);
    ASSERT_FALSE(uint32_t(~u32));
    bit_utils::make_range_mask(u64, 0, 64);
    ASSERT_FALSE(uint64_t(~u64));

    /*
     * make all range masks and check every bit in the resulting integer
     */
    size_t nbit;

    nbit = 8;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit_utils::make_range_mask(u8, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit_utils::get(u8, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit_utils::get(u8, ibit)); }
                else { ASSERT_TRUE(bit_utils::get(u8, ibit)); }
            }
        }
    }

    nbit = 16;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit_utils::make_range_mask(u16, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit_utils::get(u16, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit_utils::get(u16, ibit)); }
                else { ASSERT_TRUE(bit_utils::get(u16, ibit)); }
            }
        }
    }

    nbit = 32;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit_utils::make_range_mask(u32, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit_utils::get(u32, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit_utils::get(u32, ibit)); }
                else { ASSERT_TRUE(bit_utils::get(u32, ibit)); }
            }
        }
    }

    nbit = 64;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit_utils::make_range_mask(u64, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit_utils::get(u64, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit_utils::get(u64, ibit)); }
                else { ASSERT_TRUE(bit_utils::get(u64, ibit)); }
            }
        }
    }
}

TEST(Utils, Exsigs) {
    using namespace exsig_utils;
    /*
     * assert that all excitation signatures are encoded and decoded correctly
     */
    size_t exsig;
    for (size_t ncref = 0ul; ncref <= defs::exsig_nop_mask_frm; ++ncref) {
        for (size_t nannf = 0ul; nannf <= defs::exsig_nop_mask_frm; ++nannf) {
            for (size_t ncreb = 0ul; ncreb <= defs::exsig_nop_mask_bos; ++ncreb) {
                for (size_t nannb = 0ul; nannb <= defs::exsig_nop_mask_bos; ++nannb) {
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