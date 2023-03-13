//
// Created by Robert J. Anderson on 05/08/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/Reduction.h"

TEST(Reduction, NdAllSum) {
    const uint_t nrow = 5, ncol = 6;
    reduction::NdArray<uint_t, 2> reduction({nrow, ncol});
    for (uint_t irow = 0ul; irow < nrow; ++irow)
        for (uint_t icol = 0ul; icol < ncol; ++icol)
            reduction.m_local[{irow, icol}] = hash::in_range({irow, icol, mpi::irank()}, 0, 120);

    reduction.all_sum();
    for (uint_t irow = 0ul; irow < nrow; ++irow) {
        for (uint_t icol = 0ul; icol < ncol; ++icol) {
            uint_t tot = 0ul;
            for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank)
                tot += hash::in_range({irow, icol, irank}, 0, 120);
            auto chk = reduction.m_reduced[{irow, icol}];
            ASSERT_EQ(tot, chk);
        }
    }
}


TEST(Reduction, NdAllSumSyndicate) {
    const uint_t nrow = 5, ncol = 6;
    reduction::NdArray<uint_t, 2> reduction({nrow, ncol});
    reduction::NdArray<uint_t, 2> reduction2({nrow, ncol});
    for (uint_t irow = 0ul; irow < nrow; ++irow) {
        for (uint_t icol = 0ul; icol < ncol; ++icol) {
            reduction.m_local[{irow, icol}] = hash::in_range({irow, icol, mpi::irank()}, 0, 120);
            reduction2.m_local[{irow, icol}] = hash::in_range({irow, icol, mpi::irank()}, 0, 150);
        }
    }

    auto chk_sum_all_fn = [&](const reduction::NdArray<uint_t, 2>& r, hash::digest_t hi) {
        for (uint_t irow = 0ul; irow < nrow; ++irow) {
            for (uint_t icol = 0ul; icol < ncol; ++icol) {
                uint_t tot = 0ul;
                for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank)
                    tot += hash::in_range({irow, icol, irank}, 0, hi);
                auto chk = r.m_reduced[{irow, icol}];
                ASSERT_EQ(tot, chk);
            }
        }
    };

    v_t<reduction::Base*> syndicate = {&reduction, &reduction2};
    reduction::all_sum(syndicate);
    chk_sum_all_fn(reduction, 120);
    chk_sum_all_fn(reduction2, 150);
    reduction::all_sum(syndicate);
    chk_sum_all_fn(reduction, 120);
    chk_sum_all_fn(reduction2, 150);
    reduction::zero_local(syndicate);
    ASSERT_TRUE(reduction.m_local.is_zero());
    ASSERT_TRUE(reduction2.m_local.is_zero());

}

TEST(Reduction, CyclicScalar){
    reduction::cyclic::Scalar<uint_t, false> v;
    ASSERT_EQ(v.delta(), 0ul);
    ASSERT_EQ(v.total(), 0ul);
//
//    auto local = hash::in_range(mpi::irank(), 13, 34);
//    v.delta() = local;
//    v.all_sum();
//    int total_reduced = 0;
//    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) total_reduced+=hash::in_range(irank, 13, 34);
////    ASSERT_EQ(v.total(), local);
//    ASSERT_EQ(v.total().m_reduced, total_reduced);
}