//
// Created by Robert J. Anderson on 13/07/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/SharedArray.h"
#include "M7_lib/util/Hash.h"

TEST(SharedArray, Scalar) {
    const uint_t nelement_per_rank = 10000;
    const uint_t nelement = nelement_per_rank*mpi::nrank();
    SharedArray<double> array(nelement);

    if (mpi::on_node_i_am_root()) array.set(nelement-1, 1234);
    mpi::barrier_on_node();
    ASSERT_EQ(array[nelement-1], 1234ul);
}

TEST(SharedArray, VectorTest) {
    const uint_t nrow = 10;
    const uint_t nelement_per_rank = 10;
    v_t<SharedArray<hash::digest_t>> arrays{};
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
}
