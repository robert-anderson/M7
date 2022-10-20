//
// Created by Robert J. Anderson on 14/04/2021.
//

#include <M7_lib/connection/Connections.h>
#include "M7_lib/bilinear/FermionPromoter.h"
#include "gtest/gtest.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/util/Exsig.h"

TEST(FermionPromoter, Promoter1BodyDiagonal) {
    const uint_t nsite = 5;
    const uint_t nop_insert = 1;
    buffered::FrmOnv src(nsite);
    buffered::FrmOnv dst(nsite);
    FrmOps com(nsite);

    src = {1, 3, 4, 6, 7, 9};
    dst = {1, 3, 4, 6, 7, 9};

    conn::FrmOnv conn(src);
    conn.connect(src, dst, com);

    ASSERT_FALSE(conn.phase(src));
    FermionPromoter fp(com.size(), conn.exsig(), nop_insert);

    const auto exsig = conn.ranksig(nop_insert);
    ASSERT_EQ(exsig, exsig::ex_1100);
    buffered::MaeInds inds(exsig);
    /*
     * all diagonal promotion phases should be false (i.e. no fermi phase change)
     */
    for (uint_t icomb=0ul; icomb<fp.m_ncomb; ++icomb) {
        auto phase = fp.apply(0, conn, com, inds.m_frm);
        ASSERT_FALSE(phase);
        ASSERT_EQ(inds.m_frm.m_cre[0], com[0]);
        ASSERT_EQ(inds.m_frm.m_ann[0], com[0]);
    }
}

TEST(FermionPromoter, Promoter2BodyDiagonal) {
    const uint_t nsite = 5;
    const uint_t nop_insert = 2;
    buffered::FrmOnv src(nsite);
    buffered::FrmOnv dst(nsite);
    FrmOps com(nsite);

    src = {1, 3, 4, 6, 7, 9};
    dst = {1, 3, 4, 6, 7, 9};

    conn::FrmOnv conn(src);
    conn.connect(src, dst, com);

    ASSERT_FALSE(conn.phase(src));
    FermionPromoter fp(com.size(), conn.exsig(), nop_insert);

    const auto exsig = conn.ranksig(nop_insert);
    ASSERT_EQ(exsig, exsig::ex_2200);
    buffered::MaeInds inds(exsig);
    /*
     * all diagonal promotion phases should be false (i.e. no fermi phase change)
     */
    using namespace basic_foreach::rtnd;
    Ordered<> foreach_comb(com.size(), nop_insert);
    uint_t icomb = 0ul;
    auto fn = [&](const inds_t& insert_inds) {
        // enumerates all ways to choose nop_insert elements of the common array
        auto phase = fp.apply(icomb, conn, com, inds.m_frm);
        ASSERT_FALSE(phase);
        // since this is a diagonal element, all indices are "inserted" and should be found src order src the com array
        for (uint_t iop = 0ul; iop < inds.m_frm.m_cre.size(); ++iop)
            ASSERT_EQ(inds.m_frm.m_cre[iop], com[insert_inds[iop]]);
        for (uint_t iop = 0ul; iop < inds.m_frm.m_ann.size(); ++iop)
            ASSERT_EQ(inds.m_frm.m_ann[iop], com[insert_inds[iop]]);
        ++icomb;
    };
    foreach_comb.loop(fn);
}


TEST(FermionPromoter, Promoter2BodySingle) {
    const uint_t nsite = 5;
    const uint_t nop_insert = 1;
    buffered::FrmOnv src(nsite);
    buffered::FrmOnv dst(nsite);
    FrmOps com(nsite);

    src = {1, 3, 4, 6, 7, 9};
    dst = {1, 4, 6, 7, 8, 9};

    conn::FrmOnv conn(src);
    conn.connect(src, dst, com);

    // 3 -> 8 excitation moves through 3 occupied SQ ops => fermi phase true
    ASSERT_TRUE(conn.phase(src));
    ASSERT_EQ(conn.m_ann[0], 3);
    ASSERT_EQ(conn.m_cre[0], 8);
    FermionPromoter fp(com.size(), conn.exsig(), nop_insert);

    const auto exsig = conn.ranksig(nop_insert);
    ASSERT_EQ(exsig, exsig::ex_double);
    buffered::MaeInds inds(exsig);

    // common: 1 4 6 7 9
    bool phase;
    phase = fp.apply(0, conn, com, inds.m_frm);
    ASSERT_FALSE(phase);
    ASSERT_EQ(inds.m_frm.m_ann[0], 1);
    ASSERT_EQ(inds.m_frm.m_ann[1], 3);
    ASSERT_EQ(inds.m_frm.m_cre[0], 1);
    ASSERT_EQ(inds.m_frm.m_cre[1], 8);

    phase = fp.apply(1, conn, com, inds.m_frm);
    ASSERT_TRUE(phase);
    ASSERT_EQ(inds.m_frm.m_ann[0], 3);
    ASSERT_EQ(inds.m_frm.m_ann[1], 4);
    ASSERT_EQ(inds.m_frm.m_cre[0], 4);
    ASSERT_EQ(inds.m_frm.m_cre[1], 8);

    phase = fp.apply(2, conn, com, inds.m_frm);
    ASSERT_TRUE(phase);
    ASSERT_EQ(inds.m_frm.m_ann[0], 3);
    ASSERT_EQ(inds.m_frm.m_ann[1], 6);
    ASSERT_EQ(inds.m_frm.m_cre[0], 6);
    ASSERT_EQ(inds.m_frm.m_cre[1], 8);

    phase = fp.apply(3, conn, com, inds.m_frm);
    ASSERT_TRUE(phase);
    ASSERT_EQ(inds.m_frm.m_ann[0], 3);
    ASSERT_EQ(inds.m_frm.m_ann[1], 7);
    ASSERT_EQ(inds.m_frm.m_cre[0], 7);
    ASSERT_EQ(inds.m_frm.m_cre[1], 8);

    phase = fp.apply(4, conn, com, inds.m_frm);
    ASSERT_FALSE(phase);
    ASSERT_EQ(inds.m_frm.m_ann[0], 3);
    ASSERT_EQ(inds.m_frm.m_ann[1], 9);
    ASSERT_EQ(inds.m_frm.m_cre[0], 8);
    ASSERT_EQ(inds.m_frm.m_cre[1], 9);
}

TEST(FermionPromoter, Promoter2BodyDouble) {
    /*
     * simple test of the edge-case where no promotion is actually performed, but we need to ensure that the connection's
     * contents are faithfully reproduced src the single contributing key (fields::FermionMevInds object) emitted.
     */
    const uint_t nsite = 5;
    const uint_t nop_insert = 0;
    buffered::FrmOnv src(nsite);
    buffered::FrmOnv dst(nsite);
    FrmOps com(nsite);

    src = {1, 3, 4, 6, 7, 9};
    dst = {1, 2, 4, 6, 8, 9};

    conn::FrmOnv conn(src);
    conn.connect(src, dst, com);

    // 3, 7 -> 2, 8 excitation moves through 0 occupied SQ ops => fermi phase false
    ASSERT_FALSE(conn.phase(src));
    ASSERT_EQ(conn.m_ann[0], 3);
    ASSERT_EQ(conn.m_ann[1], 7);
    ASSERT_EQ(conn.m_cre[0], 2);
    ASSERT_EQ(conn.m_cre[1], 8);
    FermionPromoter fp(com.size(), conn.exsig(), nop_insert);

    const auto exsig = conn.ranksig(nop_insert);
    ASSERT_EQ(exsig, exsig::ex_2200);
    ASSERT_EQ(exsig, conn.exsig());
    buffered::MaeInds inds(exsig);

    bool phase = fp.apply(0, conn, com, inds.m_frm);
    ASSERT_FALSE(phase);
    ASSERT_EQ(inds.m_frm.m_ann[0], 3);
    ASSERT_EQ(inds.m_frm.m_ann[1], 7);
    ASSERT_EQ(inds.m_frm.m_cre[0], 2);
    ASSERT_EQ(inds.m_frm.m_cre[1], 8);

    ASSERT_EQ(fp.m_ncomb, 1ul);
}


TEST(FermionPromoter, Promoter3BodySingle) {
    const uint_t nsite = 5;
    const uint_t nop_insert = 2;
    buffered::FrmOnv src(nsite);
    buffered::FrmOnv dst(nsite);
    FrmOps com(nsite);

    src = {1, 4, 5, 6, 7, 9};
    dst = {1, 4, 6, 7, 8, 9};

    conn::FrmOnv conn(src);
    conn.connect(src, dst, com);

    // 5 -> 8 excitation moves through 2 occupied SQ ops => fermi phase false
    ASSERT_FALSE(conn.phase(src));
    ASSERT_EQ(conn.m_ann[0], 5);
    ASSERT_EQ(conn.m_cre[0], 8);
    FermionPromoter fp(com.size(), conn.exsig(), nop_insert);

    // number of pairs in common
    ASSERT_EQ(fp.m_ncomb, 10ul);

    const auto exsig = conn.ranksig(nop_insert);
    ASSERT_EQ(exsig, exsig::ex_triple);
    buffered::MaeInds inds(exsig);

    // common: 1 4 6 7 9
    v_t<std::pair<bool, uintv_t>> correct = {
        {0, {1, 4, 5,  1, 4, 8}},
        {1, {1, 5, 6,  1, 6, 8}},
        {1, {4, 5, 6,  4, 6, 8}},
        {1, {1, 5, 7,  1, 7, 8}},
        {1, {4, 5, 7,  4, 7, 8}},
        {0, {5, 6, 7,  6, 7, 8}},
        {0, {1, 5, 9,  1, 8, 9}},
        {0, {4, 5, 9,  4, 8, 9}},
        {1, {5, 6, 9,  6, 8, 9}},
        {1, {5, 7, 9,  7, 8, 9}},
    };

    bool phase;
    for (uint_t icomb = 0ul; icomb < correct.size(); ++icomb) {
        phase = fp.apply(icomb, conn, com, inds.m_frm);
        const auto& bench_phase = correct[icomb].first;
        const auto& bench_inds = correct[icomb].second;
        ASSERT_EQ(inds.m_frm.m_ann[0], bench_inds[0]);
        ASSERT_EQ(inds.m_frm.m_ann[1], bench_inds[1]);
        ASSERT_EQ(inds.m_frm.m_ann[2], bench_inds[2]);
        ASSERT_EQ(inds.m_frm.m_cre[0], bench_inds[3]);
        ASSERT_EQ(inds.m_frm.m_cre[1], bench_inds[4]);
        ASSERT_EQ(inds.m_frm.m_cre[2], bench_inds[5]);
        ASSERT_EQ(bench_phase, phase);
    }
}