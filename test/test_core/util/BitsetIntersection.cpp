//
// Created by Robert John Anderson on 16/07/2023.
//

#include "test_core/defs.h"
#include "M7_lib/util/BitsetIntersection.h"


TEST(BitsetIntersection, TwoSetIntersection) {
    const uint_t nelem = 1000;
    const uint_t sparsity = 3;
    const auto nbit = nelem * sparsity;
    const auto v1 = hash::unique_in_range<uint_t>(0, nelem, 0, nbit, true);
    const auto v2 = hash::unique_in_range<uint_t>(1, nelem, 0, nbit, true);
    uintv_t v3;
    std::set_intersection(v1.cbegin(), v1.cend(), v2.cbegin(), v2.cend(), std::inserter(v3, v3.begin()));
    auto b1 = bitset_isect::make_bitset(v1, nbit);
    auto b2 = bitset_isect::make_bitset(v2, nbit);
    auto b_isect = bitset_isect::isect(b1, b2);
    const auto b3 = bitset_isect::to_std_vector(b_isect);
    ASSERT_EQ(v3, b3);
}


TEST(BitsetIntersection, EnumerateAllUniqueIsects) {
    const uint_t nset = 15;
    const uint_t nelem = 200;
    const uint_t sparsity = 10;
    const uint_t nbit = nelem * sparsity;

    v_t<uintv_t> vs;
    v_t<uintv_t> bitsets;
    v_t<v_t<uintp_t>> isects;
    for (uint_t iset = 0ul; iset < nset; ++iset) {
        vs.emplace_back(hash::unique_in_range<uint_t>(iset, nelem, 0, nbit, true));
        bitsets.emplace_back(bitset_isect::make_bitset(vs.back(), nbit));
        isects.emplace_back(bitset_isect::bitset_to_siv(bitsets.back()));
    }

    auto fn = [&](const uintv_t& isets, const v_t<uintp_t>& isect) {
        uintv_t chk = vs[isets[0]];
        for (auto it = isets.cbegin() + 1; it != isets.cend(); ++it) {
            uintv_t tmp = {};
            std::set_intersection(chk.cbegin(), chk.cend(), vs[*it].cbegin(), vs[*it].cend(), std::inserter(tmp, tmp.begin()));
            chk = tmp;
        }
        ASSERT_EQ(chk, bitset_isect::to_std_vector(isect));
    };

    bitset_isect::foreach_unique(bitsets, isects, 5, true, fn);
}



