//
// Created by RJA on 24/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/basis/PermanentConnection.h"
#include "src/core/basis/Permanent.h"

TEST(PermanentConnection, Occupation){
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    // permanent -> nmode-length array of integers <= occ_cutoff
    Permanent ket(nmode, occ_cutoff);
    Permanent bra(nmode, occ_cutoff);

    ket(0) = 9;
    ket(1) = 4;
    ket(2) = 5;
    ket(3) = 12;

    bra(0) = 9;
    bra(1) = 4;
    bra(2) = 5;
    bra(3) = 13;

    //for (auto& item : ket.to_vector()) std::cout << (int)item << std::endl;

    // TODO setup occs of ket and bra
    PermanentConnection pc(ket, bra);
    ASSERT_EQ(pc.nchanged_mode(), 1);
    ASSERT_EQ(pc.changed_modes()[0], 3);
    ASSERT_EQ(pc.changes()[pc.changed_modes()[0]], -1);

    bra(3) = 10;
    pc.connect(ket, bra);
    ASSERT_EQ(pc.changes()[pc.changed_modes()[0]], 2);

    // TODO things to test... pc(k,b) == pc(b,k)
}