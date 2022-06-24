//
// Created by rja on 12/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Bit.h"

TEST(UtilBit, SetBits) {
    defs::uintv_t setbits;
    std::string bitstr;

    bitstr = "10011110"; // 158
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint8_t u8 = 0;
    for (const auto& bit: setbits) bit::set(u8, bit);
    ASSERT_EQ(u8, 158);

    bitstr = "11110111000100"; // 15812
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint16_t u16 = 0;
    for (const auto& bit: setbits) bit::set(u16, bit);
    ASSERT_EQ(u16, 15812);

    bitstr = "10100010111000001110110011101"; // 341581213
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint32_t u32 = 0;
    for (const auto& bit: setbits) bit::set(u32, bit);
    ASSERT_EQ(u32, 341581213);

    bitstr = "10000001111110111101111011110000000001100111110101001100011110"; // 2341581243132367646
    setbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='1') setbits.push_back(i);
    uint64_t u64 = 0;
    for (const auto& bit: setbits) bit::set(u64, bit);
    ASSERT_EQ(u64, 2341581243132367646);
}


TEST(UtilBit, ClrBits) {
    defs::uintv_t clrbits;
    std::string bitstr;

    bitstr = "10011110"; // 158
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint8_t u8 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit::set(u8, ibit);
    for (const auto& bit: clrbits) bit::clr(u8, bit);
    ASSERT_EQ(u8, 158);

    bitstr = "11110111000100"; // 15812
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint16_t u16 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit::set(u16, ibit);
    for (const auto& bit: clrbits) bit::clr(u16, bit);
    ASSERT_EQ(u16, 15812);

    bitstr = "10100010111000001110110011101"; // 341581213
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint32_t u32 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit::set(u32, ibit);
    for (const auto& bit: clrbits) bit::clr(u32, bit);
    ASSERT_EQ(u32, 341581213);

    bitstr = "10000001111110111101111011110000000001100111110101001100011110"; // 2341581243132367646
    clrbits.clear();
    for (size_t i=0ul; i<bitstr.size(); ++i) if (bitstr[bitstr.size()-1-i]=='0') clrbits.push_back(i);
    uint64_t u64 = 0;
    // set all bits first
    for (size_t ibit=0; ibit<bitstr.size(); ++ibit) bit::set(u64, ibit);
    for (const auto& bit: clrbits) bit::clr(u64, bit);
    ASSERT_EQ(u64, 2341581243132367646);
}

TEST(UtilBit, BitRangeMasks) {
    uint8_t u8 = 0;
    uint16_t u16 = 0;
    uint32_t u32 = 0;
    uint64_t u64 = 0;
    /*
     * whole integer edge cases: negations should be zero
     */
    bit::make_range_mask(u8, 0, 8);
    ASSERT_FALSE(uint8_t(~u8));
    bit::make_range_mask(u16, 0, 16);
    ASSERT_FALSE(uint16_t(~u16));
    bit::make_range_mask(u32, 0, 32);
    ASSERT_FALSE(uint32_t(~u32));
    bit::make_range_mask(u64, 0, 64);
    ASSERT_FALSE(uint64_t(~u64));

    /*
     * make all range masks and check every bit in the resulting integer
     */
    size_t nbit;

    nbit = 8;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit::make_range_mask(u8, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit::get(u8, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit::get(u8, ibit)); }
                else { ASSERT_TRUE(bit::get(u8, ibit)); }
            }
        }
    }

    nbit = 16;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit::make_range_mask(u16, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit::get(u16, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit::get(u16, ibit)); }
                else { ASSERT_TRUE(bit::get(u16, ibit)); }
            }
        }
    }

    nbit = 32;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit::make_range_mask(u32, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit::get(u32, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit::get(u32, ibit)); }
                else { ASSERT_TRUE(bit::get(u32, ibit)); }
            }
        }
    }

    nbit = 64;
    for (size_t ibegin = 0ul; ibegin<nbit; ++ibegin){
        for (size_t iend = ibegin+1; iend<=nbit; ++iend) {
            bit::make_range_mask(u64, ibegin, iend);
            for (size_t ibit = 0ul; ibit<nbit; ++ibit) {
                if (ibit<ibegin) { ASSERT_FALSE(bit::get(u64, ibit)); }
                else if (ibit>=iend) { ASSERT_FALSE(bit::get(u64, ibit)); }
                else { ASSERT_TRUE(bit::get(u64, ibit)); }
            }
        }
    }
}