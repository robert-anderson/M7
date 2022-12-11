//
// Created by Robert J. Anderson on 31/03/2022.
//

#include "gtest/gtest.h"

#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/connection/Connections.h"

TEST(FrmBosOnvConnection, DetectExsig) {
    const sys::frm::Basis frm_basis(6);
    const sys::bos::Basis bos_basis(6);
    buffered::FrmBosOnv src(frm_basis, bos_basis);
    src.m_frm = {{0, 1, 2}, {0, 1, 2}};
    auto dst = src;
    dst.m_bos[5] = 1;
    conn::FrmBosOnv conn(src);
    ASSERT_EQ(conn.exsig(), opsig::c_zero);
    conn.connect(src, dst);
    ASSERT_EQ(conn.exsig(), opsig::c_0010);
    conn.connect(dst, src);
    ASSERT_EQ(conn.exsig(), opsig::c_0001);
}
