//
// Created by Robert J. Anderson on 05/08/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/parallel/Reduction.h"

TEST(Reduction, NdAllSum) {
    const size_t nrow = 5, ncol = 6;
    NdReduction<size_t, 2> reduction({nrow, ncol});
    for (size_t irow = 0ul; irow < nrow; ++irow)
        for (size_t icol = 0ul; icol < ncol; ++icol)
            reduction.m_local[{irow, icol}] = utils::hash::in_range({irow, icol, mpi::irank()}, 0, 120);

    reduction.all_sum();
    for (size_t irow = 0ul; irow < nrow; ++irow) {
        for (size_t icol = 0ul; icol < ncol; ++icol) {
            size_t tot = 0ul;
            for (size_t irank = 0ul; irank < mpi::nrank(); ++irank)
                tot += utils::hash::in_range({irow, icol, irank}, 0, 120);
            auto chk = reduction.m_reduced[{irow, icol}];
            ASSERT_EQ(tot, chk);
        }
    }
}


TEST(Reduction, NdAllSumSyndicate) {
    const size_t nrow = 5, ncol = 6;
    NdReduction<size_t, 2> reduction({nrow, ncol});
    for (size_t irow = 0ul; irow < nrow; ++irow)
        for (size_t icol = 0ul; icol < ncol; ++icol)
            reduction.m_local[{irow, icol}] = utils::hash::in_range({irow, icol, mpi::irank()}, 0, 120);
    NdReduction<int, 2> reduction2({nrow, ncol});
    for (size_t irow = 0ul; irow < nrow; ++irow)
        for (size_t icol = 0ul; icol < ncol; ++icol)
            reduction2.m_local[{irow, icol}] = utils::hash::in_range({irow, icol, mpi::irank()}, 0, 150);

    ReductionSyndicate syndicate;
    syndicate.add_members(reduction, reduction2);
    syndicate.all_sum();

    for (size_t irow = 0ul; irow < nrow; ++irow) {
        for (size_t icol = 0ul; icol < ncol; ++icol) {
            size_t tot = 0ul;
            for (size_t irank = 0ul; irank < mpi::nrank(); ++irank)
                tot += utils::hash::in_range({irow, icol, irank}, 0, 120);
            auto chk = reduction.m_reduced[{irow, icol}];
            ASSERT_EQ(tot, chk);
        }
    }
    for (size_t irow = 0ul; irow < nrow; ++irow) {
        for (size_t icol = 0ul; icol < ncol; ++icol) {
            size_t tot = 0ul;
            for (size_t irank = 0ul; irank < mpi::nrank(); ++irank)
                tot += utils::hash::in_range({irow, icol, irank}, 0, 150);
            auto chk = reduction2.m_reduced[{irow, icol}];
            ASSERT_EQ(tot, chk);
        }
    }
}

