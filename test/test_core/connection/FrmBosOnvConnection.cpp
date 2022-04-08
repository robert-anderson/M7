//
// Created by rja on 31/03/2022.
//

#include "gtest/gtest.h"

#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/connection/Connections.h"

TEST(FrmBosOnvConnection, DetectExsig) {
    buffered::FrmBosOnv src(6, 6);
    src.m_frm = {{0, 1, 2}, {0, 1, 2}};
    auto dst = src;
    dst.m_bos[5] = 1;
    conn::FrmBosOnv conn(src);
    ASSERT_EQ(conn.exsig(), 0ul);
    conn.connect(src, dst);
    ASSERT_EQ(conn.exsig(), exsig_utils::ex_0010);
    conn.connect(dst, src);
    ASSERT_EQ(conn.exsig(), exsig_utils::ex_0001);
}
