//
// Created by Robert J. Anderson on 13/07/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/SharedArray.h"
#include "M7_lib/util/Hash.h"

TEST(SharedArray, VectorTest) {
    const uint_t nrow = 10;
    const uint_t nelement_per_rank = 10;
    v_t<SharedArray<hash::digest_t>> arrays;
    for (uint_t i = 0ul; i < nrow; ++i) {
        arrays.emplace_back(mpi::nrank() * nelement_per_rank);
        auto& array = arrays.back();
        for (uint_t ielem=0ul; ielem < array.size(); ++ielem) {
            array.set(ielem, hash::in_range({i, ielem}, 2, 13));
        }
    }
    for (uint_t i = 0ul; i < nrow; ++i) {
        const auto& array = arrays[i];
        for (uint_t ielem=0ul; ielem < array.size(); ++ielem) {
            ASSERT_EQ(array[ielem], hash::in_range({i, ielem}, 2, 13));
        }
    }
    auto arrays_cpy = arrays;
    for (uint_t i = 0ul; i < nrow; ++i) {
        const auto& array = arrays_cpy[i];
        for (uint_t ielem=0ul; ielem < array.size(); ++ielem) {
            ASSERT_EQ(array[ielem], hash::in_range({i, ielem}, 2, 13));
        }
    }
    auto arrays_mov = std::move(arrays);
    for (uint_t i = 0ul; i < nrow; ++i) {
        const auto& array = arrays_mov[i];
        for (uint_t ielem=0ul; ielem < array.size(); ++ielem) {
            ASSERT_EQ(array[ielem], hash::in_range({i, ielem}, 2, 13));
        }
    }
}
