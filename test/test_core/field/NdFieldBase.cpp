//
// Created by Robert J. Anderson on 06/04/2021.
//

#include <M7_lib/table/BufferedTable.h>
#include "gtest/gtest.h"
#include "M7_lib/field/Fields.h"


struct TestRow : Row {
    field::FrmOnv m_det;
    TestRow() : m_det(this, 7){}
};

TEST(NdFieldBase, Test) {

    TestRow r;
    ASSERT_EQ(r.m_size, defs::nbyte_word);

    BufferedTable<TestRow> t("", {{}});
    ASSERT_EQ(t.row_size(), defs::nbyte_word);
    t.push_back(2);
    auto r1 = t.m_row;
    r1.restart();
    auto r2 = t.m_row;
    r2.restart();

    r1.m_det = {1, 3};

    std::cout << t.to_string() << std::endl;

}