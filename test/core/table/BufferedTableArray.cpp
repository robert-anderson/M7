//
// Created by rja on 11/11/2020.
//

#include <src/core/field/Row.h>
#include <src/core/field/Fields.h>
#include <src/core/table/BufferedFields.h>
#include "gtest/gtest.h"
#include "src/core/table/BufferedTableArray.h"

TEST(BufferedTableArray, Resize){
    struct TestRow : Row {
        field::FrmOnv m_onv;
        TestRow(size_t nsite): m_onv(this, nsite){}
    };
    BufferedTableArray<TestRow> bta("Test table", 3, {{6}});
    bta.expand(4);

    bta[0].push_back(3);
    bta[0].m_row.restart();
    bta[0].m_row.m_onv = {1, 2, 3};
    bta[0].m_row.step();
    bta[0].m_row.m_onv = {1, 2, 4};
    bta[0].m_row.step();
    bta[0].m_row.m_onv = {1, 2, 5};

    bta[1].push_back(3);
    bta[1].m_row.restart();
    bta[1].m_row.m_onv = {6, 2, 3};
    bta[1].m_row.step();
    bta[1].m_row.m_onv = {6, 2, 4};
    bta[1].m_row.step();
    bta[1].m_row.m_onv = {6, 2, 5};

    bta[2].push_back(3);
    bta[2].m_row.restart();
    bta[2].m_row.m_onv = {7, 2, 3};
    bta[2].m_row.step();
    bta[2].m_row.m_onv = {7, 2, 4};
    bta[2].m_row.step();
    bta[2].m_row.m_onv = {7, 2, 5};


    buffered::FrmOnv onv(6);
    ASSERT_TRUE(onv.m_row);
    ASSERT_EQ(onv.m_table.m_hwm, 1);

    bta[0].m_row.restart();
    onv = {1, 2, 3}; ASSERT_EQ(bta[0].m_row.m_onv, onv);
    bta[0].m_row.step();
    onv = {1, 2, 4}; ASSERT_EQ(bta[0].m_row.m_onv, onv);
    bta[0].m_row.step();
    onv = {1, 2, 5}; ASSERT_EQ(bta[0].m_row.m_onv, onv);

    bta[1].m_row.restart();
    onv = {6, 2, 3}; ASSERT_EQ(bta[1].m_row.m_onv, onv);
    bta[1].m_row.step();
    onv = {6, 2, 4}; ASSERT_EQ(bta[1].m_row.m_onv, onv);
    bta[1].m_row.step();
    onv = {6, 2, 5}; ASSERT_EQ(bta[1].m_row.m_onv, onv);

    bta[2].m_row.restart();
    onv = {7, 2, 3}; ASSERT_EQ(bta[2].m_row.m_onv, onv);
    bta[2].m_row.step();
    onv = {7, 2, 4}; ASSERT_EQ(bta[2].m_row.m_onv, onv);
    bta[2].m_row.step();
    onv = {7, 2, 5}; ASSERT_EQ(bta[2].m_row.m_onv, onv);


    bta.expand(5);
    bta[0].m_row.restart();
    onv = {1, 2, 3}; ASSERT_EQ(bta[0].m_row.m_onv, onv);
    bta[0].m_row.step();
    onv = {1, 2, 4}; ASSERT_EQ(bta[0].m_row.m_onv, onv);
    bta[0].m_row.step();
    onv = {1, 2, 5}; ASSERT_EQ(bta[0].m_row.m_onv, onv);

    bta[1].m_row.restart();
    onv = {6, 2, 3}; ASSERT_EQ(bta[1].m_row.m_onv, onv);
    bta[1].m_row.step();
    onv = {6, 2, 4}; ASSERT_EQ(bta[1].m_row.m_onv, onv);
    bta[1].m_row.step();
    onv = {6, 2, 5}; ASSERT_EQ(bta[1].m_row.m_onv, onv);

    bta[2].m_row.restart();
    onv = {7, 2, 3}; ASSERT_EQ(bta[2].m_row.m_onv, onv);
    bta[2].m_row.step();
    onv = {7, 2, 4}; ASSERT_EQ(bta[2].m_row.m_onv, onv);
    bta[2].m_row.step();
    onv = {7, 2, 5}; ASSERT_EQ(bta[2].m_row.m_onv, onv);

}