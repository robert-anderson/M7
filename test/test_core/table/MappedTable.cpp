//
// Created by Robert J. Anderson on 19/07/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedTable.h"
#include "M7_lib/table/BufferedFields.h"

namespace mapped_table_test {
    struct KeyOnlyRow : Row {
        field::Number<uint_t> m_key;

        KeyOnlyRow() : m_key(this) {}

        field::Number<uint_t> &key_field() {
            return m_key;
        }
    };

    typedef buffered::MappedTable<KeyOnlyRow> key_only_table_t;

    struct ExtendedRow : KeyOnlyRow {
        field::String m_text;
        field::Numbers<float, 1> m_numbers;

        ExtendedRow() :
            m_text(this, 20, "some text"),
            m_numbers(this, {6}, "some numbers") {}
    };

    typedef buffered::MappedTable<ExtendedRow> extended_table_t;

    uint_t nitem_in_bucket(const MappedTableBase &table, uint_t ibucket) {
        const auto &bucket = table.m_buckets[ibucket];
        return std::distance(bucket.begin(), bucket.end());
    }
}

TEST(MappedTable, Empty) {
    using namespace mapped_table_test;
    key_only_table_t table("test", {});
    buffered::Number<uint_t> key;
    key = 100;
    ASSERT_FALSE(table.lookup(key));
}

TEST(MappedTable, Remap) {
    using namespace mapped_table_test;
    MappedTableOptions mapping_opts;
    mapping_opts.m_nbucket_init = 3ul;
    mapping_opts.m_remap_ratio = 0.5;
    mapping_opts.m_remap_nlookup = 10ul;
    key_only_table_t table("test", {}, mapping_opts);
    ASSERT_EQ(table.nbucket(), mapping_opts.m_nbucket_init);

    table.set_expansion_factor(1);
    buffered::Number<uint_t> key;
    key = 100;
    while ((key++) < 120) table.insert(key);

    ASSERT_EQ(table.nbucket(), mapping_opts.m_nbucket_init);
    ASSERT_EQ(nitem_in_bucket(table, 0), 6);
    ASSERT_EQ(nitem_in_bucket(table, 1), 6);
    ASSERT_EQ(nitem_in_bucket(table, 2), 8);
    ASSERT_FALSE(table.remap_due());

    ASSERT_EQ(table.m_nlookup_total, 0);
    ASSERT_EQ(table.m_nskip_total, 0);
    /*
     * bucket contents:
     * 117  114  113  111  104  102
     * 119  116  112  110  107  101
     * 120  118  115  109  108  106  105  103
     */
    key = 120;
    table.lookup(key);
    ASSERT_EQ(table.m_nlookup_total, 1);
    // 120 was the last item to be added, so no skips were required to access it
    ASSERT_EQ(table.m_nskip_total, 0);

    key = 118;
    table.lookup(key);
    ASSERT_EQ(table.m_nlookup_total, 2);
    // 118 is item 6 in bucket 2, so it should have taken 8-6-1 = 1 skip to access it
    ASSERT_EQ(table.m_nskip_total, 1);

    key = 106;
    table.lookup(key);
    ASSERT_EQ(table.m_nlookup_total, 3);
    // 106 is item 2 in bucket 2, so it should have taken 8-2-1 = 5 skips to access it
    // and we have 1 skip already so 6 in total
    ASSERT_EQ(table.m_nskip_total, 6);

    key = 110;
    table.lookup(key);
    ASSERT_EQ(table.m_nlookup_total, 4);
    // 110 is item 2 in bucket 1, so it should have taken 6-2-1 = 3 skips to access it
    // and we have 6 skips already so 9 in total
    ASSERT_EQ(table.m_nskip_total, 9);

    // remap is not yet due since there have only been 4 accesses
    ASSERT_FALSE(table.remap_due());
    key = 120;
    while (table.m_nlookup_total < table.m_mapping_opts.m_remap_nlookup) table.lookup(key);
    ASSERT_EQ(table.m_nlookup_total, table.m_mapping_opts.m_remap_nlookup);
    // remap is now due since there have been enough total lookups and skips/lookups ratio exceeds thresh
    ASSERT_TRUE(table.remap_due());
    // do the remap
    auto ratio = table.skip_lookup_ratio();
    table.attempt_remap();

    // these counters should have been reset
    ASSERT_FALSE(table.m_nlookup_total);
    ASSERT_FALSE(table.m_nskip_total);

    // remap can't be due since there haven't been any lookups since the last remap
    ASSERT_FALSE(table.remap_due());

    // make the required number of (non-skipping) lookups for remapping to be due
    key = 120;
    for (uint_t i=0ul; i<table.m_mapping_opts.m_remap_nlookup; ++i) table.lookup(key);
    // remap still shouldn't be due since there were no skips
    ASSERT_FALSE(table.remap_due());

    const auto nitem = table.nrecord();
    ASSERT_EQ(nitem, 20);
    const uint_t nbucket = mapping_opts.m_nbucket_init *
            (ratio / table.m_mapping_opts.m_remap_ratio) * (1.0 + table.get_expansion_factor());
    ASSERT_EQ(nbucket, table.nbucket());

    // check that all elements are still mapped and present after remap operation
    key = 100;
    while ((key++) < 120) ASSERT_TRUE(table.lookup(key));
}

TEST(MappedTable, Copy) {
    using namespace mapped_table_test;
    extended_table_t table("test", {});
    table.set_expansion_factor(1);
    const uint_t ibegin = 16;
    const uint_t nstep = 7;
    const uint_t nentry = 10;
    const strv_t strings {
            "Lorem ipsum", "dolor sit amet,", "consectetur", "adipiscing elit,", "sed do",
            "eiusmod tempor", "incididunt", "ut", "labore et", "dolore", "magna aliqua."
    };
    const v_t<float> floats = {1.3, 1.4, 2.5, 2.6, 3.7, 3.8};

    buffered::Number<uint_t> key;
    for (uint_t i=0ul; i<nentry; ++i) {
        key = ibegin+i*nstep;
        auto& row = table.insert(key);
        row.m_text = strings[i];
        row.m_numbers = floats;
        row.m_numbers += float(i);
    }

    /*
     * copy ctor
     */
    auto cpy = table;
    ASSERT_EQ(cpy, table);

    extended_table_t other("other", {});
    /*
     * copy assigment
     */
    other = table;
    ASSERT_EQ(other, cpy);
}