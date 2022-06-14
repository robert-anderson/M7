//
// Created by Robert J. Anderson on 19/07/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedTable.h"
#include "M7_lib/table/BufferedFields.h"

namespace mapped_table_test {
    struct MyRow : Row {
        field::Number<size_t> m_key;

        MyRow() : m_key(this) {}

        field::Number<size_t> &key_field() {
            return m_key;
        }
    };

    typedef BufferedTable<MyRow, true> table_t;

    size_t nitem_in_bucket(const MappedTableBase &table, size_t ibucket) {
        const auto &bucket = table.m_buckets[ibucket];
        return std::distance(bucket.begin(), bucket.end());
    }
}

TEST(MappedTable, Empty) {
    using namespace mapped_table_test;
    table_t table("test", {{}});
    buffered::Number<size_t> key;
    key = 100;
    auto lookup = table[key];
    ASSERT_FALSE(lookup);
}

TEST(MappedTable, Remap) {
    using namespace mapped_table_test;
    const size_t nbucket_init = 3ul;
    const size_t nlookup_remap = 10ul;
    const double remap_ratio = 0.5;
    table_t table("test", {{}, nbucket_init, nlookup_remap, remap_ratio});
    table.set_expansion_factor(1);
    ASSERT_EQ(remap_ratio, table.m_remap_ratio);
    buffered::Number<size_t> key;
    key = 100;
    while ((key++) < 120) table.insert(key);

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
    table[key];
    ASSERT_EQ(table.m_nlookup_total, 1);
    // 120 was the last item to be added, so no skips were required to access it
    ASSERT_EQ(table.m_nskip_total, 0);

    key = 118;
    table[key];
    ASSERT_EQ(table.m_nlookup_total, 2);
    // 118 is item 6 in bucket 2, so it should have taken 8-6-1 = 1 skip to access it
    ASSERT_EQ(table.m_nskip_total, 1);

    key = 106;
    table[key];
    ASSERT_EQ(table.m_nlookup_total, 3);
    // 106 is item 2 in bucket 2, so it should have taken 8-2-1 = 5 skips to access it
    // and we have 1 skip already so 6 in total
    ASSERT_EQ(table.m_nskip_total, 6);

    key = 110;
    table[key];
    ASSERT_EQ(table.m_nlookup_total, 4);
    // 110 is item 2 in bucket 1, so it should have taken 6-2-1 = 3 skips to access it
    // and we have 6 skips already so 9 in total
    ASSERT_EQ(table.m_nskip_total, 9);

    // remap is not yet due since there have only been 4 accesses
    ASSERT_FALSE(table.remap_due());
    key = 120;
    while (table.m_nlookup_total < nlookup_remap) table[key];
    ASSERT_EQ(table.m_nlookup_total, nlookup_remap);
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
    for (size_t i=0ul; i<nlookup_remap; ++i) table[key];
    // remap still shouldn't be due since there were no skips
    ASSERT_FALSE(table.remap_due());

    const auto nitem = table.nrow_nonzero();
    ASSERT_EQ(nitem, 20);
    const size_t nbucket = nbucket_init * (ratio / remap_ratio) * (1.0 + table.get_expansion_factor());
    ASSERT_EQ(nbucket, table.nbucket());

    // check that all elements are still mapped and present after remap operation
    key = 100;
    while ((key++) < 120) ASSERT_TRUE(table[key]);
}