//
// Created by Robert J. Anderson on 28/03/2022.
//

#include "gtest/gtest.h"

#include "M7_lib/basis/Suites.h"

TEST(Suites, Copy){
    const sys::frm::Basis frm_basis(6, 7, 2);
    const sys::bos::Basis bos_basis(4, 8);
    suite::Mbfs mbfs({frm_basis, bos_basis});
    auto copy_mbfs = mbfs;
    ASSERT_NE(mbfs.m_row.m_table, copy_mbfs.m_row.m_table);
    ASSERT_EQ(copy_mbfs.m_row.m_frm.m_basis, frm_basis);
    ASSERT_EQ(copy_mbfs.m_row.m_bos.m_basis, bos_basis);
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_frm.m_basis, frm_basis);
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_bos.m_basis, bos_basis);
    ASSERT_TRUE(copy_mbfs.m_row.in_range());
    defs::uintv_t frm_inds = {1, 3, 5};
    copy_mbfs.m_row.m_frm = frm_inds;
    ASSERT_EQ(copy_mbfs.m_row.m_frm, frm_inds);
    defs::uintv_t bos_inds = {4, 1, 2};
    copy_mbfs.m_row.m_bos = bos_inds;
    ASSERT_EQ(copy_mbfs.m_row.m_bos, bos_inds);
    copy_mbfs.m_row.m_frmbos = {frm_inds, bos_inds};
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_frm, frm_inds);
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_bos, bos_inds);
}
