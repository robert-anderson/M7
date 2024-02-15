//
// Created by Robert J. Anderson on 19/07/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/table/OpenAddressedTable.h"

namespace open_addressed_table_test {
    struct KeyOnlyRow : Row {
        field::Number<uint_t> m_key;

        KeyOnlyRow() : m_key(this) {}

        field::Number<uint_t> &key_field() {
            return m_key;
        }
    };

    typedef buffered::OpenAddressedTable<KeyOnlyRow> key_only_table_t;

    struct ExtendedRow : KeyOnlyRow {
        field::String m_text;
        field::Numbers<float, 1> m_numbers;

        ExtendedRow() :
            m_text(this, 20, "some text"),
            m_numbers(this, {6}, "some numbers") {}
    };

    typedef buffered::OpenAddressedTable<ExtendedRow> extended_table_t;
}

TEST(OpenAddressedTable, Empty) {
    using namespace open_addressed_table_test;
    key_only_table_t table("test", {}, 0.5, Owner::local());
    buffered::Number<uint_t> key;
    key = 100;
    ASSERT_FALSE(table.lookup(key));
}

TEST(OpenAddressedTable, InsertLocal) {
    using namespace open_addressed_table_test;
    key_only_table_t table("test", {}, 0.5, Owner::local());
    table.set_expansion_factor(0.5);
    buffered::Number<uint_t> key;
    auto nums = hash::unique_in_range(0, 100, 23, 340);
    // resize to 1/4 of the final number of keys so that OpenAddresser has to remap on the fly
    table.resize(nums.size()/4);
    for (const auto& num: nums) {
        key = num;
        ASSERT_FALSE(table.lookup(key));
        table.insert(key);
        ASSERT_TRUE(table.lookup(key));
    }

    for (const auto& num: nums) {
        key = num;
        ASSERT_TRUE(table.lookup(key));
    }

    // try to insert one of the same keys again
    key = nums[0];
    ASSERT_TRUE(table.lookup(key));
    table.lookup_or_insert(key);
    ASSERT_TRUE(table.lookup(key));
    ASSERT_EQ(table.nrow_in_use(), nums.size());
}

TEST(OpenAddressedTable, InsertShared) {
    using namespace open_addressed_table_test;
    key_only_table_t table("test", {}, 0.5, Owner::shared());
    table.set_expansion_factor(0.5);
    buffered::Number<uint_t> key;
    auto nums = hash::unique_in_range(0, 100, 23, 340);
    // resize to the final number of keys since resizes are synchronizing operations
    table.resize(nums.size());
    if (table.i_can_modify()) {
        for (const auto &num: nums) {
            key = num;
            ASSERT_FALSE(table.lookup(key));
            table.insert(key);
            ASSERT_TRUE(table.lookup(key));
        }
    }
    auto more_nums = hash::unique_in_range(1, 100, 340, 500);
    /*
     * resize to the final number of keys since resizes are synchronizing operations, this triggers a remap of the
     * OpenAddresser
     */
    table.resize(nums.size() + more_nums.size());
    if (table.i_can_modify()) {
        for (const auto &num: more_nums) {
            key = num;
            ASSERT_FALSE(table.lookup(key));
            table.insert(key);
            ASSERT_TRUE(table.lookup(key));
        }
    }
    table.end_sync();
    ASSERT_EQ(table.nrow_in_use(), nums.size() + more_nums.size());

    for (const auto& num: nums) {
        key = num;
        ASSERT_TRUE(table.lookup(key));
    }
    for (const auto& num: more_nums) {
        key = num;
        ASSERT_TRUE(table.lookup(key));
    }

    // try to insert one of the same keys again
    key = nums[0];
    ASSERT_TRUE(table.lookup(key));
    if (table.i_can_modify()) table.lookup_or_insert(key);
    table.end_sync();
    ASSERT_TRUE(table.lookup(key));
    ASSERT_EQ(table.nrow_in_use(), nums.size() + more_nums.size());
}