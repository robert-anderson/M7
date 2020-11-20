//
// Created by rja on 11/11/2020.
//

#include <src/core/field/Elements.h>
#include "gtest/gtest.h"
#include "src/core/table/BufferedTableArray.h"

struct TestTable : TableX {
    fields::FermionOnv m_fonv;
    TestTable(size_t nsite):
            m_fonv(this, nsite, "configuration"){}
};

TEST(BufferedTableArray, Resize){
    BufferedTableArray<TestTable> bta("Test table", 3, 6);
    bta.expand(4);

    bta[0].push_back(3);
    bta[0].m_fonv(0) = {1, 2, 3};
    bta[0].m_fonv(1) = {1, 2, 4};
    bta[0].m_fonv(2) = {1, 2, 5};

    bta[1].push_back(3);
    bta[1].m_fonv(0) = {6, 2, 3};
    bta[1].m_fonv(1) = {6, 2, 4};
    bta[1].m_fonv(2) = {6, 2, 5};

    bta[2].push_back(3);
    bta[2].m_fonv(0) = {7, 2, 3};
    bta[2].m_fonv(1) = {7, 2, 4};
    bta[2].m_fonv(2) = {7, 2, 5};

    elements::FermionOnv fonv(6);
    fonv = {1, 2, 3}; ASSERT_EQ(bta[0].m_fonv(0), fonv);
    fonv = {1, 2, 4}; ASSERT_EQ(bta[0].m_fonv(1), fonv);
    fonv = {1, 2, 5}; ASSERT_EQ(bta[0].m_fonv(2), fonv);

    fonv = {6, 2, 3}; ASSERT_EQ(bta[1].m_fonv(0), fonv);
    fonv = {6, 2, 4}; ASSERT_EQ(bta[1].m_fonv(1), fonv);
    fonv = {6, 2, 5}; ASSERT_EQ(bta[1].m_fonv(2), fonv);

    fonv = {7, 2, 3}; ASSERT_EQ(bta[2].m_fonv(0), fonv);
    fonv = {7, 2, 4}; ASSERT_EQ(bta[2].m_fonv(1), fonv);
    fonv = {7, 2, 5}; ASSERT_EQ(bta[2].m_fonv(2), fonv);


    bta.expand(5);
    fonv = {1, 2, 3}; ASSERT_EQ(bta[0].m_fonv(0), fonv);
    fonv = {1, 2, 4}; ASSERT_EQ(bta[0].m_fonv(1), fonv);
    fonv = {1, 2, 5}; ASSERT_EQ(bta[0].m_fonv(2), fonv);

    fonv = {6, 2, 3}; ASSERT_EQ(bta[1].m_fonv(0), fonv);
    fonv = {6, 2, 4}; ASSERT_EQ(bta[1].m_fonv(1), fonv);
    fonv = {6, 2, 5}; ASSERT_EQ(bta[1].m_fonv(2), fonv);

    fonv = {7, 2, 3}; ASSERT_EQ(bta[2].m_fonv(0), fonv);
    fonv = {7, 2, 4}; ASSERT_EQ(bta[2].m_fonv(1), fonv);
    fonv = {7, 2, 5}; ASSERT_EQ(bta[2].m_fonv(2), fonv);
}