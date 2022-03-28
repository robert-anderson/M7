//
// Created by rja on 28/03/2022.
//

#include "gtest/gtest.h"

#include "M7_lib/basis/Suites.h"

TEST(Suites, Copy){
    const size_t nsite = 6;
    const size_t nmode = 3;
    suite::Mbfs mbfs({nsite, nmode});
    auto copy_mbfs = mbfs;
    ASSERT_NE(mbfs.m_row.m_table, copy_mbfs.m_row.m_table);
    ASSERT_EQ(copy_mbfs.m_row.m_frm.m_nsite, nsite);
    ASSERT_EQ(copy_mbfs.m_row.m_bos.m_nmode, nmode);
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_frm.m_nsite, nsite);
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_bos.m_nmode, nmode);
    ASSERT_TRUE(copy_mbfs.m_row.in_range());
    defs::inds frm_inds = {1, 3, 5};
    copy_mbfs.m_row.m_frm = frm_inds;
    ASSERT_EQ(copy_mbfs.m_row.m_frm, frm_inds);
    defs::inds bos_inds = {4, 1, 2};
    copy_mbfs.m_row.m_bos = bos_inds;
    ASSERT_EQ(copy_mbfs.m_row.m_bos, bos_inds);
    copy_mbfs.m_row.m_frmbos = {frm_inds, bos_inds};
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_frm, frm_inds);
    ASSERT_EQ(copy_mbfs.m_row.m_frmbos.m_bos, bos_inds);
}