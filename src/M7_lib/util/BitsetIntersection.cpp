//
// Created by rja on 16/07/23.
//

#include "BitsetIntersection.h"


typedef bitset_isect::siv_t siv_t;
typedef bitset_isect::vsiv_t vsiv_t;

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

void bitset_isect::isect(const uintv_t& bitset1, const uintv_t& bitset2, uint_t iword_begin, uint_t iword_end, siv_t& siv) {
    DEBUG_ASSERT_EQ(bitset1.size(), bitset2.size(), "incompatible bitsets");
    DEBUG_ASSERT_LE(iword_begin, iword_end, "nonsensical word range");
    DEBUG_ASSERT_LT(iword_begin, bitset1.size(), "nonsensical word range");
    DEBUG_ASSERT_LE(iword_end, bitset1.size(), "nonsensical word range");
    siv.clear();
    for (uint_t iword = iword_begin; iword < iword_end; ++iword) {
        const auto isect_word = bitset1[iword] & bitset2[iword];
        if (isect_word) siv.emplace_back(iword, isect_word);
    }
}

void bitset_isect::isect(const uintv_t& bitset1, const uintv_t& bitset2, siv_t& siv) {
    isect(bitset1, bitset2, 0, bitset1.size(), siv);
}

siv_t bitset_isect::isect(const uintv_t& bitset1, const uintv_t& bitset2) {
    siv_t siv;
    isect(bitset1, bitset2, 0, bitset1.size(), siv);
    return siv;
}

siv_t bitset_isect::bitset_to_siv(const uintv_t& bitset) {
    return isect(bitset, bitset);
}

void bitset_isect::bitset_to_siv(const uintv_t& bitset, siv_t& siv) {
    isect(bitset, bitset, siv);
}

vsiv_t bitset_isect::isect_many(const v_t<uintv_t>& bitsets1, const v_t<uintv_t>& bitsets2,
                          uint_t iword_begin, uint_t iword_end) {
    vsiv_t isects;
    DEBUG_ASSERT_EQ(bitsets1.size(), bitsets2.size(), "numbers of bitsets should match");
    isects.reserve(bitsets1.size());
    while (isects.size() != isects.capacity()) {
        isects.push_back({});
        isect(bitsets1[isects.size()-1], bitsets2[isects.size()-1], iword_begin, iword_end, isects.back());
    }
    return isects;
}

vsiv_t bitset_isect::bitset_to_siv_many(const v_t<uintv_t> &bitsets) {
    return bitset_to_siv_many(bitsets, 0, bitsets[0].size());
}

vsiv_t bitset_isect::bitset_to_siv_many(const v_t<uintv_t>& bitsets, uint_t iword_begin, uint_t iword_end) {
    return isect_many(bitsets, bitsets, iword_begin, iword_end);
}

void bitset_isect::isect(siv_t& siv, const uintv_t& bitset) {
    auto ifill = 0ul;
    for (auto& pair: siv) {
        const auto isect_word = bitset[pair.first] & pair.second;
        if (isect_word) siv[ifill++] = {pair.first, isect_word};
    }
    siv.resize(ifill);
}

void bitset_isect::isect(const siv_t& siv_in, const uintv_t& bitset, siv_t& siv_out) {
    siv_out.clear();
    for (auto& pair: siv_in) {
        uint_t isect_word = pair.first < bitset.size() ? bitset[pair.first] : 0ul;
        isect_word &= pair.second;
        if (isect_word) siv_out.emplace_back(pair.first, isect_word);
    }
}

std::set<uint_t> bitset_isect::to_std_set(const siv_t& siv) {
    std::set<uint_t> s;
    foreach_in_isect(siv, [&s](uint_t i){s.insert(i);});
    return s;
}

uintv_t bitset_isect::to_std_vector(const siv_t& siv) {
    uintv_t v;
    foreach_in_isect(siv, [&v](uint_t i){v.push_back(i);});
    return v;
}
