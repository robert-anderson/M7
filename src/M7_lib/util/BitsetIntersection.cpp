//
// Created by rja on 16/07/23.
//

#include "BitsetIntersection.h"

uintv_t bitset_isect::make_bitset(const uintv_t& isetbits, uint_t nbit) {
    uintv_t bitset;
    bitset.resize(integer::divceil(nbit, sizeof(uint_t) * CHAR_BIT));
    for (auto& ind: isetbits) {
        DEBUG_ASSERT_LT(ind, nbit, "bit index OoB")
        const auto iword = ind / (sizeof(uint_t) * CHAR_BIT);
        const auto ibit = ind % (sizeof(uint_t) * CHAR_BIT);
        bit::set(bitset[iword], ibit);
    }
    return bitset;
}

uintv_t bitset_isect::make_bitset(const uintv_t& isetbits) {
    if (isetbits.empty()) return {};
    const auto it = std::max_element(isetbits.cbegin(), isetbits.cend());
    return make_bitset(isetbits, *it);
}

void bitset_isect::make_isect(const uintv_t& bitset1, const uintv_t& bitset2, uint_t iword_begin, uint_t iword_end,
                              v_t<uintp_t>& isect) {
    DEBUG_ASSERT_EQ(bitset1.size(), bitset2.size(), "incompatible bitsets");
    DEBUG_ASSERT_LE(iword_begin, iword_end, "nonsensical word range");
    DEBUG_ASSERT_LT(iword_begin, bitset1.size(), "nonsensical word range");
    DEBUG_ASSERT_LE(iword_end, bitset1.size(), "nonsensical word range");
    isect.clear();
    for (uint_t iword = iword_begin; iword < iword_end; ++iword) {
        const auto isect_word = bitset1[iword] & bitset2[iword];
        if (isect_word) isect.emplace_back(iword, isect_word);
    }
}

void bitset_isect::make_isect(const uintv_t& bitset1, const uintv_t& bitset2, v_t<uintp_t>& isect) {
    make_isect(bitset1, bitset2, 0, bitset1.size(), isect);
}

v_t<uintp_t> bitset_isect::make_isect(const uintv_t& bitset1, const uintv_t& bitset2) {
    v_t<uintp_t> isect;
    make_isect(bitset1, bitset2, 0, bitset1.size(), isect);
    return isect;
}

v_t<uintp_t> bitset_isect::make_isect(const uintv_t& bitset) {
    return make_isect(bitset, bitset);
}

void bitset_isect::make_isect(const uintv_t& bitset, v_t<uintp_t>& isect) {
    make_isect(bitset, bitset, isect);
}

v_t<v_t<uintp_t>> bitset_isect::make_isects(const v_t<uintv_t>& bitsets1, const v_t<uintv_t>& bitsets2,
                          uint_t iword_begin, uint_t iword_end) {
    v_t<v_t<uintp_t>> isects;
    DEBUG_ASSERT_EQ(bitsets1.size(), bitsets2.size(), "numbers of bitsets should match");
    isects.reserve(bitsets1.size());
    while (isects.size() != isects.capacity()) {
        isects.push_back({});
        make_isect(bitsets1[isects.size()-1], bitsets2[isects.size()-1], iword_begin, iword_end, isects.back());
    }
    return isects;
}

v_t<v_t<uintp_t>> bitset_isect::make_isects(const v_t<uintv_t>& bitsets, uint_t iword_begin, uint_t iword_end) {
    return make_isects(bitsets, bitsets, iword_begin, iword_end);
}

void bitset_isect::isect(v_t<uintp_t>& isect, const uintv_t& bitset) {
    auto ifill = 0ul;
    for (auto& pair: isect) {
        const auto isect_word = bitset[pair.first] & pair.second;
        if (isect_word) isect[ifill++] = {pair.first, isect_word};
    }
    isect.resize(ifill);
}

void bitset_isect::make_isect(const v_t<uintp_t>& isect_in, const uintv_t& bitset, v_t<uintp_t>& isect_out) {
    isect_out.clear();
    for (auto& pair: isect_in) {
        uint_t isect_word = pair.first < bitset.size() ? bitset[pair.first] : 0ul;
        isect_word &= pair.second;
        if (isect_word) isect_out.emplace_back(pair.first, isect_word);
    }
}

std::set<uint_t> bitset_isect::to_std_set(const v_t<uintp_t>& isect) {
    std::set<uint_t> s;
    foreach_in_isect(isect, [&s](uint_t i){s.insert(i);});
    return s;
}

uintv_t bitset_isect::to_std_vector(const v_t<uintp_t>& isect) {
    uintv_t v;
    foreach_in_isect(isect, [&v](uint_t i){v.push_back(i);});
    return v;
}