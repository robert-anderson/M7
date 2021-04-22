//
// Created by rja on 14/04/2021.
//

#include <src/core/basis/Connections.h>
#include <src/core/observables/MevGroup.h>
#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"

TEST(MevGroup, Promoter) {

    const size_t nsite = 5;
    const size_t nop_insert = 2;
    buffered::FermionOnv in(nsite);
    buffered::FermionOnv out(nsite);

    in = {1, 3, 4, 5, 6, 7};
    out = {0, 1, 3, 5, 6, 7};

    conn::Antisym<0> conn(nsite);

    conn.connect(in, out);

    ASSERT_FALSE(conn.phase());
    ASSERT_EQ(conn.ann(0), 4);
    ASSERT_EQ(conn.cre(0), 0);

    FermionPromoter fp(conn.ncom(), nop_insert);
    for (size_t icomb=0ul; icomb<fp.m_ncomb; ++icomb){
        for (size_t i=0ul; i<nop_insert; ++i) {
            std::cout << (size_t) fp.m_all_combs[icomb * nop_insert + i] << " ";
        }
        std::cout <<std::endl;
    }

    buffered::FermionMevInds inds(conn.nexcit()+nop_insert);
    ASSERT_EQ(inds.m_row->m_fields[1]->m_row_offset, 3ul);

    auto phase = fp.apply(2, conn, inds);
    std::cout << inds << std::endl;
    std::cout << phase << std::endl;


}

/*
TEST(MevGroup, Test){
    const size_t nsite = 6;
    buffered::FermionOnv bra(nsite);
    buffered::FermionOnv ket(nsite);

    bra = {1, 4, 5, 7, 8, 11};
    std::cout << bra.to_string() << std::endl;
    const size_t iann = 8;
    const size_t icre = 6;
    ket = bra;
    ket.clr(iann);
    ket.set(icre);
    conn::Antisym<> conn(nsite);

    conn.connect(bra, ket);
    ASSERT_EQ(conn.nexcit(), 1);
    ASSERT_EQ(conn.ann(0), iann);
    ASSERT_EQ(conn.cre(0), icre);

    conn.zero();
    //conn.add(iann, 1, icre, 1);
    conn.add(1, iann, 1, icre);
    conn.apply(bra, ket);

    ASSERT_TRUE(conn.phase());
}
 */