//
// Created by rja on 14/04/2021.
//

#include <src/core/connection/Connections.h>
#include <src/core/observables/MevGroup.h>
#include <src/core/util/Foreach.h>
#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"

TEST(MevGroup, Promoter1BodyDiagonal) {
    const size_t nsite = 5;
    const size_t nop_insert = 1;
    buffered::FrmOnv in(nsite);
    buffered::FrmOnv out(nsite);
    FrmOps com(nsite);

    in = {1, 3, 4, 6, 7, 9};
    out = {1, 3, 4, 6, 7, 9};

    conn::FrmOnv conn(nsite);
    conn.connect(in, out, com);

    ASSERT_FALSE(conn.phase(in));
    FermionPromoter fp(com.size(), nop_insert);

    buffered::FermionMevInds inds(conn.m_cre.size()+nop_insert);
    /*
     * all diagonal promotion phases should be false (i.e. no fermi phase change)
     */
    for (size_t icomb=0ul; icomb<fp.m_ncomb; ++icomb) {
        auto phase = fp.apply(0, conn, com, inds);
        ASSERT_FALSE(phase);
        ASSERT_EQ(inds.m_cre[0], com[0]);
        ASSERT_EQ(inds.m_ann[0], com[0]);
    }
}

TEST(MevGroup, Promoter2BodyDiagonal) {
    const size_t nsite = 5;
    const size_t nop_insert = 2;
    buffered::FrmOnv in(nsite);
    buffered::FrmOnv out(nsite);
    FrmOps com(nsite);

    in = {1, 3, 4, 6, 7, 9};
    out = {1, 3, 4, 6, 7, 9};

    conn::FrmOnv conn(nsite);
    conn.connect(in, out, com);

    ASSERT_FALSE(conn.phase(in));
    FermionPromoter fp(com.size(), nop_insert);

    const auto rank = conn.m_cre.size()+nop_insert;
    ASSERT_EQ(rank, 2);
    buffered::FermionMevInds inds(rank);
    /*
     * all diagonal promotion phases should be false (i.e. no fermi phase change)
     */
    foreach::rtnd::Ordered<> foreach_comb(com.size(), nop_insert);
    size_t icomb = 0ul;
    auto fn = [&]() {
        auto phase = fp.apply(icomb, conn, com, inds);
        ASSERT_FALSE(phase);
        for (size_t iop = 0ul; iop < rank; ++iop) {
            ASSERT_EQ(inds.m_ann[iop], com[foreach_comb[iop]]);
            ASSERT_EQ(inds.m_cre[iop], com[foreach_comb[iop]]);
        }
        ++icomb;
    };
    foreach_comb(fn);
}


TEST(MevGroup, Promoter2BodySingle) {
    const size_t nsite = 5;
    const size_t nop_insert = 1;
    buffered::FrmOnv in(nsite);
    buffered::FrmOnv out(nsite);
    FrmOps com(nsite);

    in = {1, 3, 4, 6, 7, 9};
    out = {1, 4, 6, 7, 8, 9};

    conn::FrmOnv conn(nsite);
    conn.connect(in, out, com);

    // 3 -> 8 excitation moves through 3 occupied SQ op inds => fermi phase -1
    ASSERT_TRUE(conn.phase(in));
    ASSERT_EQ(conn.m_ann[0], 3);
    ASSERT_EQ(conn.m_cre[0], 8);
    FermionPromoter fp(com.size(), nop_insert);

    const auto rank = conn.m_cre.size()+nop_insert;
    ASSERT_EQ(rank, 2);
    buffered::FermionMevInds inds(rank);

    // common: 1 4 6 7 9
    bool phase;
    phase = fp.apply(0, conn, com, inds);
    ASSERT_FALSE(phase);
    ASSERT_EQ(inds.m_ann[0], 1);
    ASSERT_EQ(inds.m_ann[1], 3);
    ASSERT_EQ(inds.m_cre[0], 1);
    ASSERT_EQ(inds.m_cre[1], 8);

    phase = fp.apply(1, conn, com, inds);
    ASSERT_TRUE(phase);
    ASSERT_EQ(inds.m_ann[0], 3);
    ASSERT_EQ(inds.m_ann[1], 4);
    ASSERT_EQ(inds.m_cre[0], 4);
    ASSERT_EQ(inds.m_cre[1], 8);

    phase = fp.apply(2, conn, com, inds);
    ASSERT_TRUE(phase);
    ASSERT_EQ(inds.m_ann[0], 3);
    ASSERT_EQ(inds.m_ann[1], 6);
    ASSERT_EQ(inds.m_cre[0], 6);
    ASSERT_EQ(inds.m_cre[1], 8);

    phase = fp.apply(3, conn, com, inds);
    ASSERT_TRUE(phase);
    ASSERT_EQ(inds.m_ann[0], 3);
    ASSERT_EQ(inds.m_ann[1], 7);
    ASSERT_EQ(inds.m_cre[0], 7);
    ASSERT_EQ(inds.m_cre[1], 8);

    phase = fp.apply(4, conn, com, inds);
    ASSERT_FALSE(phase);
    ASSERT_EQ(inds.m_ann[0], 3);
    ASSERT_EQ(inds.m_ann[1], 9);
    ASSERT_EQ(inds.m_cre[0], 8);
    ASSERT_EQ(inds.m_cre[1], 9);
}

TEST(MevGroup, Promoter2BodyDouble) {
    /*
     * simple test of the edge-case where no promotion is actually performed, but we need to ensure that the connection's
     * contents are faithfully reproduced in the single contributing key (fields::FermionMevInds object) emitted.
     */
    const size_t nsite = 5;
    const size_t nop_insert = 0;
    buffered::FrmOnv in(nsite);
    buffered::FrmOnv out(nsite);
    FrmOps com(nsite);

    in = {1, 3, 4, 6, 7, 9};
    out = {1, 2, 4, 6, 8, 9};

    conn::FrmOnv conn(nsite);
    conn.connect(in, out, com);

    // 3, 7 -> 2, 8 excitation moves through 0 occupied SQ op inds => fermi phase +1
    ASSERT_FALSE(conn.phase(in));
    ASSERT_EQ(conn.m_ann[0], 3);
    ASSERT_EQ(conn.m_ann[1], 7);
    ASSERT_EQ(conn.m_cre[0], 2);
    ASSERT_EQ(conn.m_cre[1], 8);
    FermionPromoter fp(com.size(), nop_insert);

    const auto rank = conn.m_cre.size()+nop_insert;
    ASSERT_EQ(rank, 2);
    buffered::FermionMevInds inds(rank);

    bool phase = fp.apply(0, conn, com, inds);
    ASSERT_FALSE(phase);
    ASSERT_EQ(inds.m_ann[0], 3);
    ASSERT_EQ(inds.m_ann[1], 7);
    ASSERT_EQ(inds.m_cre[0], 2);
    ASSERT_EQ(inds.m_cre[1], 8);
}


TEST(MevGroup, FermionRdm){
    //FermionRdm
}