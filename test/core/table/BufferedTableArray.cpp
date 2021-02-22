//
// Created by rja on 11/11/2020.
//

#include <src/core/field/Row.h>
#include <src/core/field/Fields.h>
#include <src/core/field/BufferedFields.h>
#include "gtest/gtest.h"
#include "src/core/field/BufferedTableArray.h"

TEST(BufferedTableArray, Resize){
    struct TestRow : Row {
        fields::Onv<0> m_fonv;
        TestRow(size_t nsite): m_fonv(this, nsite){}
    };
    BufferedTableArray<TestRow> bta("Test table", 3, 6);
    bta.expand(4);

    bta[0].push_back(3);
    bta[0].m_row.restart();
    bta[0].m_row.m_fonv = {1, 2, 3};
    bta[0].m_row.step();
    bta[0].m_row.m_fonv = {1, 2, 4};
    bta[0].m_row.step();
    bta[0].m_row.m_fonv = {1, 2, 5};

    bta[1].push_back(3);
    bta[1].m_row.restart();
    bta[1].m_row.m_fonv = {6, 2, 3};
    bta[1].m_row.step();
    bta[1].m_row.m_fonv = {6, 2, 4};
    bta[1].m_row.step();
    bta[1].m_row.m_fonv = {6, 2, 5};

    bta[2].push_back(3);
    bta[2].m_row.restart();
    bta[2].m_row.m_fonv = {7, 2, 3};
    bta[2].m_row.step();
    bta[2].m_row.m_fonv = {7, 2, 4};
    bta[2].m_row.step();
    bta[2].m_row.m_fonv = {7, 2, 5};


    buffered::FermionOnv fonv(6);
    ASSERT_TRUE(fonv.m_row);
    ASSERT_EQ(fonv.m_table.m_hwm, 1);

    bta[0].m_row.restart();
    fonv = {1, 2, 3}; ASSERT_EQ(bta[0].m_row.m_fonv, fonv);
    bta[0].m_row.step();
    fonv = {1, 2, 4}; ASSERT_EQ(bta[0].m_row.m_fonv, fonv);
    bta[0].m_row.step();
    fonv = {1, 2, 5}; ASSERT_EQ(bta[0].m_row.m_fonv, fonv);

    bta[1].m_row.restart();
    fonv = {6, 2, 3}; ASSERT_EQ(bta[1].m_row.m_fonv, fonv);
    bta[1].m_row.step();
    fonv = {6, 2, 4}; ASSERT_EQ(bta[1].m_row.m_fonv, fonv);
    bta[1].m_row.step();
    fonv = {6, 2, 5}; ASSERT_EQ(bta[1].m_row.m_fonv, fonv);

    bta[2].m_row.restart();
    fonv = {7, 2, 3}; ASSERT_EQ(bta[2].m_row.m_fonv, fonv);
    bta[2].m_row.step();
    fonv = {7, 2, 4}; ASSERT_EQ(bta[2].m_row.m_fonv, fonv);
    bta[2].m_row.step();
    fonv = {7, 2, 5}; ASSERT_EQ(bta[2].m_row.m_fonv, fonv);


    bta.expand(5);
    bta[0].m_row.restart();
    fonv = {1, 2, 3}; ASSERT_EQ(bta[0].m_row.m_fonv, fonv);
    bta[0].m_row.step();
    fonv = {1, 2, 4}; ASSERT_EQ(bta[0].m_row.m_fonv, fonv);
    bta[0].m_row.step();
    fonv = {1, 2, 5}; ASSERT_EQ(bta[0].m_row.m_fonv, fonv);

    bta[1].m_row.restart();
    fonv = {6, 2, 3}; ASSERT_EQ(bta[1].m_row.m_fonv, fonv);
    bta[1].m_row.step();
    fonv = {6, 2, 4}; ASSERT_EQ(bta[1].m_row.m_fonv, fonv);
    bta[1].m_row.step();
    fonv = {6, 2, 5}; ASSERT_EQ(bta[1].m_row.m_fonv, fonv);

    bta[2].m_row.restart();
    fonv = {7, 2, 3}; ASSERT_EQ(bta[2].m_row.m_fonv, fonv);
    bta[2].m_row.step();
    fonv = {7, 2, 4}; ASSERT_EQ(bta[2].m_row.m_fonv, fonv);
    bta[2].m_row.step();
    fonv = {7, 2, 5}; ASSERT_EQ(bta[2].m_row.m_fonv, fonv);

}