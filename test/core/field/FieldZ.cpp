//
// Created by rja on 09/02/2021.
//

#include <src/core/fieldz/MappedTableZ.h>
#include <src/core/fieldz/BufferedTableZ.h>
#include "gtest/gtest.h"
#include "src/core/fieldz/FieldsZ.h"
#include "src/core/fieldz/TableZ.h"
#include "src/core/fieldz/BufferedFields.h"

TEST(FieldZ, Copying) {
    struct TestRow : RowZ {
        fieldsz::Vector<float> m_vector;
        fieldsz::Vectors<float> m_vectors;

        TestRow() : m_vector(this, 9), m_vectors(this, 4, 5){}
    };

    TestRow row;
    ASSERT_EQ(&row.m_vector, row.m_fields[0]);
    ASSERT_EQ(row.m_vector.m_row_offset, 0ul);
    ASSERT_EQ(row.m_fields.size(), 2ul);
    TestRow rowcopy(row);
    ASSERT_FALSE(rowcopy.m_fields.empty());
    /*
     * check that all pointers have been updated correctly
     */
    ASSERT_EQ(&rowcopy.m_vector, rowcopy.m_fields[0]);
    ASSERT_EQ(rowcopy.m_vector.m_row, &rowcopy);
    ASSERT_EQ(rowcopy.m_fields.size(), 2ul);
    ASSERT_EQ(rowcopy.m_vector.m_row_offset, 0ul);
}

TEST(FieldZ, Test) {
    struct TestRow : RowZ {
        fieldsz::Vectors<float> m_array;
        fieldsz::FermiBosOnv m_fbonv;
        fieldsz::Flags m_flags;

        TestRow() :
                m_array(this, 3, 14),
                m_fbonv(this, 12),
                m_flags(this, 7) {}

        fieldsz::FermiBosOnv &m_key_field = m_fbonv;
    };

    BufferedTableZ<TestRow> table("Test", {});

    ASSERT_EQ(&table.m_row, table.m_row.m_array.m_row);
    ASSERT_EQ(table.m_row.m_array.m_row_offset, 0ul);
    ASSERT_EQ(table.m_row.m_flags.m_item_size, 4);
    ASSERT_EQ(table.m_row.m_flags.m_size, 4);

    const size_t nrow = 30;
    table.push_back(nrow);
    auto& row = table.m_row;
    size_t irow = 0ul;
    for (row.restart(); row.in_range(); row.step()) irow++;

    ASSERT_EQ(irow, nrow);
}