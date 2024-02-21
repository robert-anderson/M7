//
// Created by Robert John Anderson on 12/02/2024.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedTable.h"
#include "M7_lib/table/OpenAddresser.h"

const buf_t* get_key(const strv_t& keys, uint_t i){
    return reinterpret_cast<const buf_t*>(keys[i].data());
}

TEST(OpenAddresser, LocalInsertAndLookup) {
    const uint_t nchar = 12;
    using row_t = SingleFieldRow<StringField>;
    buffered::Table<row_t> table({nchar}, Owner::local());
    const double fmax = 0.3;
    Buffer buffer(1, Owner::local());
    OpenAddresser oa(table, 0, table.m_row.m_field.m_size, fmax);
    oa.set_buffer(&buffer);
    strv_t keys = {
        "lorem_______", "ipsum_______", "dolor_______", "sit_________", "amet,_______", "consectetur_",
        "adipiscing__", "elit,_______", "sed_________", "do__________", "eiusmod_____", "tempor______",
        "incididunt__", "ut__________", "labore______", "et__________", "dolore______", "magna_______",
        "aliqua______"
    };
    table.push_back(keys.size());
    oa.resize();

    // the row positions in table at which each key is mapped
    uintv_t current_entries(keys.size(), ~0ul);

    auto row = table.m_row;
    const uint_t n_move = 400;
    // each move is either an insertion or an erasure
    for (uint_t i_move = 0; i_move < n_move; ++i_move) {
        const auto move = hash::in_range(i_move, 0, 2 * keys.size());
        const auto i_key = move / 2;
        const auto do_insert = move % 2;
        if (do_insert) {
            // insertion
            if (current_entries[i_key] == ~0ul) {
                // not already inserted
                // scan to a free slot in the table
                for (row.restart(); row && !row.m_field.is_clear(); ++row){}
                ASSERT_TRUE(row.is_deref_valid());
                row.m_field = keys[i_key];
                ASSERT_TRUE(oa.insert(row.index()));
                current_entries[i_key] = row.index();
            }
            else {
                // already inserted
                const auto ind = oa.lookup(get_key(keys, i_key));
                ASSERT_EQ(ind, current_entries[i_key]);
            }
        }
        else {
            // erasure
            if (current_entries[i_key] != ~0ul) {
                // already inserted: can be erased
                const auto ind = oa.lookup(get_key(keys, i_key));
                ASSERT_NE(ind, ~0ul);
                ASSERT_EQ(ind, current_entries[i_key]);
                oa.erase(ind);
                row.jump(ind);
                row.m_field.clear();
                current_entries[i_key] = ~0ul;
            }
        }
    }
    for (uint_t i = 0; i < keys.size(); ++i) ASSERT_EQ(oa.lookup(get_key(keys, i)), current_entries[i]);

    str_t key = "not_real_key";
    ASSERT_EQ(oa.lookup(reinterpret_cast<const buf_t*>(key.data())), ~0ul);
}

TEST(OpenAddresser, SharedInsertAndLookup) {
    const uint_t nchar = 12;
    using row_t = SingleFieldRow<StringField>;
    buffered::Table<row_t> table({nchar}, Owner::shared(0));
    const double fmax = 0.3;
    OpenAddresser oa(table, table.owner(), 0, table.m_row.m_field.m_size, fmax);
    strv_t keys = {
            "lorem_______", "ipsum_______", "dolor_______", "sit_________", "amet,_______", "consectetur_",
            "adipiscing__", "elit,_______", "sed_________", "do__________", "eiusmod_____", "tempor______",
            "incididunt__", "ut__________", "labore______", "et__________", "dolore______", "magna_______",
            "aliqua______"
    };
    table.push_back(keys.size());
    oa.resize();

    // the row positions in table at which each key is mapped
    uintv_t current_entries(keys.size(), ~0ul);

    auto row = table.m_row;
    const uint_t n_move = 400;
    // each move is either an insertion or an erasure
    if (mpi::g_irank_shmem == 0) {
        for (uint_t i_move = 0; i_move < n_move; ++i_move) {
            const auto move = hash::in_range(i_move, 0, 2 * keys.size());
            const auto i_key = move / 2;
            const auto do_insert = move % 2;
            if (do_insert) {
                // insertion
                if (current_entries[i_key] == ~0ul) {
                    // not already inserted
                    // scan to a free slot in the table
                    for (row.restart(); row && !row.m_field.is_clear(); ++row){}
                    ASSERT_TRUE(row.is_deref_valid());
                    row.m_field = keys[i_key];
                    ASSERT_TRUE(oa.insert(row.index()));
                    current_entries[i_key] = row.index();
                }
                else {
                    // already inserted
                    const auto ind = oa.lookup(get_key(keys, i_key));
                    ASSERT_EQ(ind, current_entries[i_key]);
                }
            }
            else {
                // erasure
                if (current_entries[i_key] != ~0ul) {
                    // already inserted: can be erased
                    const auto ind = oa.lookup(get_key(keys, i_key));
                    ASSERT_NE(ind, ~0ul);
                    ASSERT_EQ(ind, current_entries[i_key]);
                    oa.erase(ind);
                    row.jump(ind);
                    row.m_field.clear();
                    current_entries[i_key] = ~0ul;
                }
            }
        }
    }
    mpi::bcast(current_entries, 0);

    for (uint_t i = 0; i < keys.size(); ++i) ASSERT_EQ(oa.lookup(get_key(keys, i)), current_entries[i]);
    mpi::barrier();
    str_t key = "not_real_key";
    ASSERT_EQ(oa.lookup(reinterpret_cast<const buf_t*>(key.data())), ~0ul);
}
