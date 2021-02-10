//
// Created by rja on 09/02/2021.
//

#include "gtest/gtest.h"
#include "src/core/fieldz/FieldsZ.h"
#include "src/core/fieldz/TableZ.h"

TEST(FieldZ, Copying) {
    struct TestRow : RowZ {
        fieldsz::NumberArrays<float, 1, 2> m_array;
        TestRow() : m_array(this, {1}, {2, 3}) {}
    };

    TestRow row;
    ASSERT_EQ(&row.m_array.get<0>(), row.m_fields[0]);
    ASSERT_EQ(row.m_array.get<0>().m_row_offset, 0ul);
    TestRow rowcopy(row);
    ASSERT_FALSE(rowcopy.m_fields.empty());
    /*
     * check that all pointers have been updated correctly
     */
    ASSERT_EQ(&rowcopy.m_array.get<0>(), rowcopy.m_fields[0]);
    ASSERT_EQ(rowcopy.m_array.get<0>().m_row, &rowcopy);
    ASSERT_EQ(rowcopy.m_array.m_row, &rowcopy);
    ASSERT_EQ(rowcopy.m_array.get<0>().m_row_offset, 0ul);
}


TEST(FieldZ, Test) {
    struct TestRow : RowZ {
        fieldsz::NumberArrays<float, 1, 2> m_array;
        TestRow() :
        m_array(this, {1}, {2,3}){}
    };

    Buffer buffer("Test", 1);
    TableZ<TestRow> table({});

    ASSERT_EQ(&table.m_row, table.m_row.m_array.m_row);
    ASSERT_EQ(table.m_row.m_array.get<0>().m_row_offset, 0ul);

    table.set_buffer(&buffer);
    const size_t nrow = 102;
    table.push_back(nrow);
    table.m_row.restart();
    for (size_t irow = 1ul; irow<nrow; ++irow)
        ASSERT_TRUE(table.m_row.try_step());
    ASSERT_FALSE(table.m_row.try_step());

    table.m_row.jump(23);
    table.m_row.m_array.get<0>()[3] = 1.231;
    table.m_row.jump(22);
    table.m_row.step();
    std::cout <<              table.m_row.m_array.get<0>().to_string()    << std::endl;
    std::cout <<              table.m_row.to_string()    << std::endl;
    /*
    std::cout <<
              table.to_string()
    << std::endl;
     */
}