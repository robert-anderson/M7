//
// Created by Robert John Anderson on 13/09/2023.
//

#include "test_core/defs.h"
#include "M7_lib/table/Smuvi.h"
#include "M7_lib/table/BufferedFields.h"

TEST(Smuvi, LocalTable) {
    typedef SingleFieldRow<field::Number<uint_t>> row_t;
    buffered::smuvi::LocalTable<row_t> table("test", row_t(), {});
    buffered::Number<uint_t> work;
    v_t<uintp_t> pairs = {
        {5, 11}, {3, 15}, {4, 19}, {3, 18}, {4, 18}, {5, 15}, {5, 14}, {3, 17}, {4, 17}
    };
    std::map<uint_t, uintv_t> combined_pairs;
    for (const auto& pair: pairs) {
        auto it = combined_pairs.find(pair.first);
        if (it==combined_pairs.end()) it = combined_pairs.insert({pair.first, {}}).first;
        it->second.push_back(pair.second);
    }

    for (const auto& pair : pairs) {
        work = pair.first;
        table.append(work, pair.second);
    }

    auto row = table.m_row;
    for (row.restart(); row; ++row) {
        const uint_t key = row.m_field;
        auto inds_ptr = table.inds(row);
        ASSERT_TRUE(inds_ptr);
        ASSERT_EQ(*inds_ptr, combined_pairs[key]);
    }
}