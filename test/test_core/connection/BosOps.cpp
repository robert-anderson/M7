//
// Created by anderson on 02/06/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/connection/BosOnvConnection.h"


TEST(BosOps, FromImodes1) {
    const size_t nmode = 6;
    BosOps ops(nmode);
    ops.set(4);
    ASSERT_EQ(ops.get()[0], 4);
    ASSERT_EQ(ops[0].m_imode, 4);
    ASSERT_EQ(ops[0].m_nop, 1);
}

TEST(BosOps, FromImodes2) {
    const size_t nmode = 6;
    BosOps ops(nmode);
    ops.set(4, 5);
    ASSERT_EQ(ops.get()[0], 4);
    ASSERT_EQ(ops.get()[1], 5);
    ASSERT_EQ(ops[0].m_imode, 4);
    ASSERT_EQ(ops[1].m_imode, 5);
    ASSERT_EQ(ops[0].m_nop, 1);
    ASSERT_EQ(ops[1].m_nop, 1);
    ops.set(4, 4);
    ASSERT_EQ(ops.get()[0], 4);
    ASSERT_EQ(ops.get()[1], 4);
    ASSERT_EQ(ops[0].m_imode, 4);
    ASSERT_EQ(ops[0].m_nop, 2);
}
TEST(BosOps, FromImodes3) {
    const size_t nmode = 6;
    BosOps ops(nmode);
    ops.set(1, 3, 5);
    ASSERT_EQ(ops.get()[0], 1);
    ASSERT_EQ(ops.get()[1], 3);
    ASSERT_EQ(ops.get()[2], 5);
    ASSERT_EQ(ops[0].m_imode, 1);
    ASSERT_EQ(ops[0].m_nop, 1);
    ASSERT_EQ(ops[1].m_imode, 3);
    ASSERT_EQ(ops[1].m_nop, 1);
    ASSERT_EQ(ops[2].m_imode, 5);
    ASSERT_EQ(ops[2].m_nop, 1);
    ops.set(3, 3, 5);
    ASSERT_EQ(ops.get()[0], 3);
    ASSERT_EQ(ops.get()[1], 3);
    ASSERT_EQ(ops.get()[2], 5);
    ASSERT_EQ(ops[0].m_imode, 3);
    ASSERT_EQ(ops[0].m_nop, 2);
    ASSERT_EQ(ops[1].m_imode, 5);
    ASSERT_EQ(ops[1].m_nop, 1);
    ops.set(3, 5, 5);
    ASSERT_EQ(ops.get()[0], 3);
    ASSERT_EQ(ops.get()[1], 5);
    ASSERT_EQ(ops.get()[2], 5);
    ASSERT_EQ(ops[0].m_imode, 3);
    ASSERT_EQ(ops[0].m_nop, 1);
    ASSERT_EQ(ops[1].m_imode, 5);
    ASSERT_EQ(ops[1].m_nop, 2);
    ops.set(5, 5, 5);
    ASSERT_EQ(ops.get()[0], 5);
    ASSERT_EQ(ops.get()[1], 5);
    ASSERT_EQ(ops.get()[2], 5);
    ASSERT_EQ(ops[0].m_imode, 5);
    ASSERT_EQ(ops[0].m_nop, 3);
}