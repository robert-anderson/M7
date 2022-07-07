//
// Created by Robert J. Anderson on 10/02/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedFields.h"

TEST(BitsetField, SetFromInds) {
    const uint_t nbit = 100;
    const uint_t nsetbit = 50;
    // deliberately using non-standard container type
    buffered::Bitset<uint16_t> field(nbit);
    auto setbits = convert::vector<uint_t>(hash::unique_in_range(0, nsetbit, 0, nbit, true));
    field = setbits;
    auto it = setbits.cbegin();
    for (uint_t ibit=0ul; ibit<nbit; ++ibit){
        if (it==setbits.cend() || ibit<*it) {
            ASSERT_FALSE(field.get(ibit));
        }
        else {
            ASSERT_TRUE(field.get(ibit));
            ++it;
        }
    }
}

TEST(BitsetField, SetRange) {
    const uint_t nbit = 100;
    // deliberately using non-standard container type
    buffered::Bitset<uint16_t> field(nbit);
    for (uint_t ibegin=0ul; ibegin<nbit; ++ibegin){
        for (uint_t iend=ibegin+1; iend<nbit; ++iend) {
            field.zero();
            field.set_range(ibegin, iend);
            for (uint_t ibit=0ul; ibit<nbit; ++ibit){
                if (ibit<ibegin) {ASSERT_FALSE(field.get(ibit));}
                else if (ibit>=iend) {ASSERT_FALSE(field.get(ibit));}
                else {ASSERT_TRUE(field.get(ibit));}
            }
        }
    }
}

TEST(BitsetField, ClrRange) {
    const uint_t nbit = 100;
    // deliberately using non-standard container type
    buffered::Bitset<uint16_t> field(nbit);
    for (uint_t ibegin=0ul; ibegin<nbit; ++ibegin){
        for (uint_t iend=ibegin+1; iend<nbit; ++iend) {
            field.set();
            field.clr_range(ibegin, iend);
            for (uint_t ibit=0ul; ibit<nbit; ++ibit){
                if (ibit<ibegin) {ASSERT_TRUE(field.get(ibit));}
                else if (ibit>=iend) {ASSERT_TRUE(field.get(ibit));}
                else {ASSERT_FALSE(field.get(ibit));}
            }
        }
    }
}
