//
// Created by Robert J. Anderson on 02/10/2020.
//

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/table//BufferedTable.h"
#include "gtest/gtest.h"

TEST(BufferedTable, Empty) {
    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> table({});
    ASSERT_EQ(table.capacity(), 0);
    ASSERT_EQ(table.nrow_in_use(), 0);
    ASSERT_EQ(table.m_bw.m_size, 0);
    ASSERT_EQ(table.cbegin(), nullptr);
    table.m_row.restart();
    ASSERT_FALSE(table.m_row.is_deref_valid());
    auto cpy = table;
    ASSERT_EQ(cpy.capacity(), 0);
    ASSERT_EQ(cpy.nrow_in_use(), 0);
    ASSERT_EQ(cpy.m_bw.m_size, 0);
    ASSERT_EQ(cpy.cbegin(), nullptr);
    cpy.m_row.restart();
    ASSERT_FALSE(cpy.m_row.is_deref_valid());
    ASSERT_NE(&cpy.m_row, &table.m_row);
    ASSERT_EQ(&table, table.m_row.m_table);
    ASSERT_EQ(&cpy, cpy.m_row.m_table);
}

TEST(BufferedTable, NodeShared) {
    typedef SingleFieldRow<field::Number<hash::digest_t>> row_t;
    buffered::Table<row_t> table({}, true);
    const uint_t nrow = 23;
    table.resize(nrow);
    for (uint_t irow=0ul; irow<nrow; ++irow) {
        table.m_row.push_back_jump();
        table.m_row.m_field = hash::in_range({irow, mpi::irank()}, 12, 190);
    }
    mpi::barrier_on_node();
    /*
     * check there is no data corruption from race conditions, and that the saved values are those from the node roots
     */
    for (table.m_row.restart(); table.m_row; ++table.m_row){
        ASSERT_EQ(table.m_row.m_field, hash::in_range({table.m_row.index(), mpi::my_node_root_irank()}, 12, 190));
    }
}

TEST(BufferedTable, AllGathervEmpty) {
    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> src_table("src", {});
    buffered::Table<row_t> dst_table("dst", {});
    dst_table.all_gatherv(src_table);
    ASSERT_EQ(dst_table.nrow_in_use(), 0ul);
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
    ASSERT_EQ(dst_table.nrow_in_use(), nrow_global);
    auto& row = dst_table.m_row;
    row.restart();
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank){
        auto nrow = get_nrow(irank);
        for (uint_t irow=0ul; irow<nrow; ++irow){
            ASSERT_EQ(row.m_field, get_value(irank, irow));
            ++row;
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
        ASSERT_EQ(dst_table.nrow_in_use(), nrow_global);
        auto &row = dst_table.m_row;
        row.restart();
        for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
            auto nrow = get_nrow(irank);
            for (uint_t irow = 0ul; irow < nrow; ++irow) {
                ASSERT_EQ(row.m_field, get_value(irank, irow));
                ++row;
            }
        }
    }
    else {
        ASSERT_EQ(dst_table.nrow_in_use(), 0ul);
    }
}