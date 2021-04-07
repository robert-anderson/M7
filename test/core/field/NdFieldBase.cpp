//
// Created by rja on 06/04/2021.
//

#include <src/core/table/BufferedTable.h>
#include "gtest/gtest.h"
#include "src/core/field/Fields.h"


struct TestRow : Row {
    fields::FermionOnv m_det;
    TestRow() : m_det(this, 7){}
};

TEST(NdFieldBase, Test) {

    TestRow r;
    ASSERT_EQ(r.m_dsize, 1);

    BufferedTable<TestRow> t("", {{}});
    ASSERT_EQ(t.m_row_dsize, 1);
    t.push_back(2);
    auto r1 = t.m_row;
    r1.restart();
    auto r2 = t.m_row;
    r2.restart();

    r1.m_det = {1, 3};

    std::cout <<
              t.to_string()
              << std::endl;

}