//
// Created by anderson on 1/25/22.
//

#include <src/core/table/BufferedTable.h>
#include "gtest/gtest.h"
#include "src/core/field/CompositeField.h"
#include "src/core/field/Fields.h"


struct CompFieldTestRow : Row {
    field::FrmOnv m_frm_onv;
    CompositeField<field::FrmOnv, field::FrmOnv> m_comp_field;
    CompFieldTestRow() :
            m_frm_onv({this, 6, "alfa"}),
            m_comp_field({this, 6, "bravo"}, {this, 8, "charlie"}){}
};

struct NestedCompFieldTestRow : Row {
    field::FrmOnv m_frm_onv;
    typedef CompositeField<field::FrmOnv, field::FrmOnv> frm_onv_pair_t;
    frm_onv_pair_t m_comp_field;
    CompositeField<frm_onv_pair_t, field::BosOnv> m_nested_field;

    NestedCompFieldTestRow() :
            m_frm_onv({this, 6, "alfa"}),
            m_comp_field({this, 6, "bravo"}, {this, 8, "charlie"}),
            m_nested_field({{this, 10, "delta"}, {this, 5, "echo"}},
                           {this, 9, "foxtrot"}){}
};

TEST(CompositeField, CopyMoveSemantics) {

    CompFieldTestRow row;
    ASSERT_EQ(row.nfield(), 3ul);
    for (size_t i=0ul; i<row.nfield(); ++i) ASSERT_EQ(row.m_fields[i]->m_row, &row);
    ASSERT_EQ(row.m_child, nullptr);

    auto row_copy = row;
    for (size_t i=0ul; i<row_copy.nfield(); ++i) ASSERT_EQ(row_copy.m_fields[i]->m_row, &row_copy);

    BufferedTable<CompFieldTestRow> bt("test table", {{}});
    std::cout << bt.m_row.m_fields.size() << std::endl;
    bt.resize(10);
    bt.push_back(4);
    std::cout << bt.to_string() << std::endl;
}

TEST(CompositeField, NestedCopyMoveSemantics) {

    //NestedCompFieldTestRow row;

    BufferedTable<NestedCompFieldTestRow> bt("test table", {{}});
    std::cout << bt.m_row.m_fields.size() << std::endl;
    bt.resize(10);
    bt.push_back(4);
    std::cout << bt.to_string() << std::endl;
}