//
// Created by Robert John Anderson on 2020-02-04.
//

#include <gtest/gtest.h>
#include "src/data/BitfieldNew.h"

TEST(Bitfield, InternalMemoryIntegrity) {
    BitfieldNew b(90);
    defs::inds setinds{0, 13, 64, 78, 89};
    b.set(setinds);
    size_t isetind{0ul};
    for (auto i{0ul}; i<b.m_nbit; ++i) {
        if (i==setinds[isetind]) {
            isetind++;
            ASSERT_TRUE(b.get(i));
        }
        else ASSERT_FALSE(b.get(i));
    }
}