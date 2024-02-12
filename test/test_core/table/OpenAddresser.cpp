//
// Created by Robert John Anderson on 12/02/2024.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedTable.h"
#include "M7_lib/table/OpenAddresser.h"


TEST(OpenAddresser, LocalInsertAndLookup) {
    using row_t = SingleFieldRow<StringField>;
    const uint_t nchar = 10;
    buffered::Table<row_t> table({nchar}, Owner::local());
    const double f_max = 0.5;
    OpenAddresser oa(table, 0, table.m_row.m_field.m_size, 0.5);
    auto row = table.m_row;
    row.push_back_jump();
    str_t key = "helloworld";
    row.m_field = key;
    // table has changes size, so sync the size of the open addresser
    oa.resize();
    oa.insert(0);
    ASSERT_EQ(oa.lookup(reinterpret_cast<const buf_t*>(key.c_str())), 0);
}