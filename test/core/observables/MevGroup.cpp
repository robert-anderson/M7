//
// Created by rja on 14/04/2021.
//

#include <src/core/basis/Connections.h>
#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"

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