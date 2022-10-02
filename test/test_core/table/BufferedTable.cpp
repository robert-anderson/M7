//
// Created by Robert J. Anderson on 02/10/2020.
//

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/table//BufferedTable.h"
#include "gtest/gtest.h"

TEST(BufferedTable, Empty) {
    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> table({});
    ASSERT_EQ(table.nrecord(), 0);
    ASSERT_EQ(table.m_hwm, 0);
    ASSERT_EQ(table.m_bw.m_size, 0);
    ASSERT_EQ(table.m_bw.m_begin, nullptr);
    table.m_row.restart();
    ASSERT_FALSE(table.m_row.in_range());
    ASSERT_FALSE(table.m_row.ptr_in_range());
    auto cpy = table;
    ASSERT_EQ(cpy.nrecord(), 0);
    ASSERT_EQ(cpy.m_hwm, 0);
    ASSERT_EQ(cpy.m_bw.m_size, 0);
    ASSERT_EQ(cpy.m_bw.m_begin, nullptr);
    cpy.m_row.restart();
    ASSERT_FALSE(cpy.m_row.in_range());
    ASSERT_FALSE(cpy.m_row.ptr_in_range());
    ASSERT_NE(&cpy.m_row, &table.m_row);
    ASSERT_EQ(&table, table.m_row.m_table);
    ASSERT_EQ(&cpy, cpy.m_row.m_table);
}

TEST(BufferedTable, AllGathervEmpty) {
    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> src_table("src", {});
    buffered::Table<row_t> dst_table("dst", {});
    dst_table.all_gatherv(src_table);
    ASSERT_EQ(dst_table.m_hwm, 0ul);
}

TEST(BufferedTable, AllGatherv) {
    /*
     * if there is more than one rank, have the second one (arbitrary choice) be empty to test the ability of the
     * gathering functionality to deal with nullptr buffer dbegins.
     */
    const uint_t irank_empty = mpi::nrank()==1 ? ~0ul: 1ul;
    auto get_nrow = [irank_empty](uint_t irank){return irank==irank_empty ? 0ul : hash::in_range(irank, 3, 10);};
    auto get_value = [](uint_t irank, uint_t irow){return int(hash::in_range(irow * (irank+1), 0, 100));};
    const uint_t nrow_local = get_nrow(mpi::irank());
    uint_t nrow_global = 0ul;
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) nrow_global+=get_nrow(irank);
    ASSERT_EQ(nrow_global,mpi::all_sum(nrow_local));

    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> src_table("src", {});
    src_table.resize(nrow_local);
    for (uint_t irow = 0ul; irow<nrow_local; ++irow){
        src_table.m_row.push_back_jump();
        src_table.m_row.m_field = get_value(mpi::irank(), irow);
    }
    buffered::Table<row_t> dst_table("dst", {});
    dst_table.all_gatherv(src_table);
    ASSERT_EQ(dst_table.m_hwm, nrow_global);
    auto& row = dst_table.m_row;
    row.restart();
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank){
        auto nrow = get_nrow(irank);
        for (uint_t irow=0ul; irow<nrow; ++irow){
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
    const uint_t irank_empty = mpi::nrank()==1 ? ~0ul: 1ul;
    auto get_nrow = [irank_empty](uint_t irank){return irank==irank_empty ? 0ul : hash::in_range(irank, 3, 10);};
    auto get_value = [](uint_t irank, uint_t irow){return int(hash::in_range(irow * (irank+1), 0, 100));};
    const uint_t nrow_local = get_nrow(mpi::irank());
    uint_t nrow_global = 0ul;
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) nrow_global+=get_nrow(irank);
    ASSERT_EQ(nrow_global,mpi::all_sum(nrow_local));

    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> src_table("src", {});
    src_table.resize(nrow_local);
    for (uint_t irow = 0ul; irow<nrow_local; ++irow){
        src_table.m_row.push_back_jump();
        src_table.m_row.m_field = get_value(mpi::irank(), irow);
    }
    buffered::Table<row_t> dst_table("dst", {});

    dst_table.gatherv(src_table);
    if (mpi::i_am_root()) {
        ASSERT_EQ(dst_table.m_hwm, nrow_global);
        auto &row = dst_table.m_row;
        row.restart();
        for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
            auto nrow = get_nrow(irank);
            for (uint_t irow = 0ul; irow < nrow; ++irow) {
                ASSERT_EQ(row.m_field, get_value(irank, irow));
                row.step();
            }
        }
    }
    else {
        ASSERT_EQ(dst_table.m_hwm, 0ul);
    }
}