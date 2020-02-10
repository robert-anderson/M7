//
// Created by Robert John Anderson on 2020-01-27.
//

#include <gtest/gtest.h>

#include "../../src/fermion/Determinant.h"

/*
 * This test set is much like that defined for CharLookups, except
 * here the multi-byte case is considered.
 */


/*
TEST(DeterminantConnection, NSetBits) {
    Determinant det(18);
    EXPECT_EQ(det.nsetbits(), 0ul);

    det.zero();
    det.set(defs::inds{0,1,4,5});
    EXPECT_EQ(det.nsetbits(), 4ul);

    det.zero();
    det.set(defs::inds{0,3,6,9,12,15});
    EXPECT_EQ(det.nsetbits(), 6ul);
    det.zero();
    EXPECT_EQ(det.nsetbits(), 0ul);
}*/