//
// Created by rja on 09/02/2021.
//

#include <src/core/fieldz/MappedTableZ.h>
#include <src/core/fieldz/BufferedTableZ.h>
#include "gtest/gtest.h"
#include "src/core/fieldz/FieldsZ.h"
#include "src/core/fieldz/TableZ.h"
#include "src/core/fieldz/ElementsZ.h"

TEST(FieldZ, Copying) {
    struct TestRow : RowZ {
        fieldsz::NumberArrays<float, 1, 2> m_array;
        fieldsz::NumberArrays<float, 1, 2> m_array2;

        TestRow() : m_array(this, {1}, {2, 3}),
                    m_array2(this, {1}, {2, 3}) {}
    };

    TestRow row;
    ASSERT_EQ(&row.m_array.get<0>(), row.m_fields[0]);
    ASSERT_EQ(row.m_array.get<0>().m_row_offset, 0ul);
    ASSERT_EQ(row.m_fields.size(), 2ul);
    TestRow rowcopy(row);
    ASSERT_FALSE(rowcopy.m_fields.empty());
    /*
     * check that all pointers have been updated correctly
     */
    ASSERT_EQ(&rowcopy.m_array.get<0>(), rowcopy.m_fields[0]);
    ASSERT_EQ(rowcopy.m_array.get<0>().m_row, &rowcopy);
    ASSERT_EQ(rowcopy.m_array.m_row, &rowcopy);
    ASSERT_EQ(rowcopy.m_fields.size(), 2ul);
    ASSERT_EQ(rowcopy.m_array.get<0>().m_row_offset, 0ul);
}

TEST(FieldZ, Test) {
    struct TestRow : RowZ {
        fieldsz::NumberArrays<float, 1, 2> m_array;
        fieldsz::FermiBosOnv m_fbonv;
        fieldsz::Flags<2> m_flags;

        TestRow() :
                m_array(this, {1}, {2, 3}),
                m_fbonv(this, 7),
                m_flags(this, {4, 7}) {}

        fieldsz::FermiBosOnv &m_key_field = m_fbonv;
    };

    BufferedTableZ<TestRow> table("Test", {});

    ASSERT_EQ(&table.m_row, table.m_row.m_array.m_row);
    ASSERT_EQ(table.m_row.m_array.get<0>().m_row_offset, 0ul);

    const size_t nrow = 30;
    table.push_back(nrow);
    table.m_row.restart();
    for (size_t irow = 1ul; irow < nrow; ++irow)
        ASSERT_TRUE(table.m_row.try_step());
    ASSERT_FALSE(table.m_row.try_step());

    ASSERT_EQ(table.m_row.m_flags().m_size, 4);

    MappedTableZ<TestRow> mt({}, 40);

//
//    struct FermionOnvRowZ {
//        struct WrappedRow : RowZ {
//            fieldsz::FermionOnv m_field;
//        };
//        FermionOnvRowZ(size_t nsite) : m_field(this, nsite) {}
//    };


    //elementsz::NumberArray<double, 1> det(NumberFieldZ<double, 1>({3}));
    elementsz::FermionOnv det(5);
    std::cout << det.to_string() << std::endl;

}