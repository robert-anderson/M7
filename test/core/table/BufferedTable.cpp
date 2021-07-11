//
// Created by rja on 02/10/2020.
//

#include <src/core/table/BufferedFields.h>
#include "src/core/table//BufferedTable.h"
#include "gtest/gtest.h"

TEST(BufferedTable, Empty) {
    typedef SingleFieldRow<fields::Number<int>> row_t;
    BufferedTable<row_t> table("", {{}});
    ASSERT_EQ(table.m_nrow, 0);
    ASSERT_EQ(table.m_hwm, 0);
    ASSERT_EQ(table.m_bw.dsize(), 0);
    ASSERT_EQ(table.m_bw.m_dbegin, nullptr);
    auto& row = table.m_row;
    row.restart();
    ASSERT_FALSE(row.in_range());
    ASSERT_FALSE(row.ptr_in_range());
}

TEST(BufferedTable, AllGathervEmpty) {
    typedef SingleFieldRow<fields::Number<int>> row_t;
    BufferedTable<row_t> src_table("src", {{}});
    BufferedTable<row_t> dst_table("dst", {{}});
    dst_table.all_gatherv(src_table);
    ASSERT_EQ(dst_table.m_hwm, 0ul);
}

TEST(BufferedTable, AllGatherv) {
    /*
     * if there is more than one rank, have the second one (arbitrary choice) be empty to test the ability of the
     * gathering functionality to deal with nullptr buffer dbegins.
     */
    const size_t irank_empty = mpi::nrank()==1 ? ~0ul: 1ul;
    auto get_nrow = [irank_empty](size_t irank){return irank==irank_empty ? 0ul : hashing::in_range(irank, 3, 10);};
    auto get_value = [](size_t irank, size_t irow){return hashing::in_range(irow * (irank+1), 0, 100);};
    const size_t nrow_local = get_nrow(mpi::irank());
    size_t nrow_global = 0ul;
    for (size_t irank=0ul; irank<mpi::nrank(); ++irank) nrow_global+=get_nrow(irank);
    ASSERT_EQ(nrow_global,mpi::all_sum(nrow_local));

    typedef SingleFieldRow<fields::Number<int>> row_t;
    BufferedTable<row_t> src_table("src", {{}});
    src_table.resize(nrow_local);
    for (size_t irow = 0ul; irow<nrow_local; ++irow){
        src_table.m_row.push_back_jump();
        src_table.m_row.m_field = get_value(mpi::irank(), irow);
    }
    BufferedTable<row_t> dst_table("dst", {{}});
    dst_table.all_gatherv(src_table);
    ASSERT_EQ(dst_table.m_hwm, nrow_global);
    auto& row = dst_table.m_row;
    row.restart();
    for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
        auto nrow = get_nrow(irank);
        for (size_t irow=0ul; irow<nrow; ++irow){
            ASSERT_EQ(row.m_field, get_value(irank, irow));
            row.step();
        }
    }
}

TEST(BufferedTable, Gatherv) {
    /*
     * if there is more than one rank, have the second one (arbitrary choice) be empty to test the ability of the
     * gathering functionality to deal with nullptr buffer dbegins.
     */
    const size_t irank_empty = mpi::nrank()==1 ? ~0ul: 1ul;
    auto get_nrow = [irank_empty](size_t irank){return irank==irank_empty ? 0ul : hashing::in_range(irank, 3, 10);};
    auto get_value = [](size_t irank, size_t irow){return hashing::in_range(irow * (irank+1), 0, 100);};
    const size_t nrow_local = get_nrow(mpi::irank());
    size_t nrow_global = 0ul;
    for (size_t irank=0ul; irank<mpi::nrank(); ++irank) nrow_global+=get_nrow(irank);
    ASSERT_EQ(nrow_global,mpi::all_sum(nrow_local));

    typedef SingleFieldRow<fields::Number<int>> row_t;
    BufferedTable<row_t> src_table("src", {{}});
    src_table.resize(nrow_local);
    for (size_t irow = 0ul; irow<nrow_local; ++irow){
        src_table.m_row.push_back_jump();
        src_table.m_row.m_field = get_value(mpi::irank(), irow);
    }
    BufferedTable<row_t> dst_table("dst", {{}});

    dst_table.gatherv(src_table);
    if (mpi::i_am_root()) {
        ASSERT_EQ(dst_table.m_hwm, nrow_global);
        auto &row = dst_table.m_row;
        row.restart();
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
            auto nrow = get_nrow(irank);
            for (size_t irow = 0ul; irow < nrow; ++irow) {
                ASSERT_EQ(row.m_field, get_value(irank, irow));
                row.step();
            }
        }
    }
    else {
        ASSERT_EQ(dst_table.m_hwm, 0ul);
    }
}