//
// Created by Robert John Anderson on 2020-01-31.
//

#include <gtest/gtest.h>

/*
#include "../../src/enumerator/CommonBitEnumerator.h"

TEST(CommonBitEnumerator, Bitfield) {
    Bitfield bitfield1(200);
    Bitfield bitfield2(200);
    defs::inds tmp{};
    tmp.assign({1, 7, 8, 32, 89, 123, 199});
    bitfield1.set(tmp);
    tmp.assign({1, 4, 7, 32, 89, 125, 199});
    bitfield2.set(tmp);
    CommonBitEnumerator<true> enumerator;
    enumerator = CommonBitEnumerator(bitfield1, bitfield2);

    tmp.assign({1, 7, 32, 89, 199});
    size_t setbit;
    while(enumerator.next(setbit)){
        ASSERT_TRUE(tmp[enumerator.counter()]==setbit);
    }
    ASSERT_EQ(enumerator.counter(), tmp.size());
    ASSERT_EQ(enumerator.counter(), bitfield1.nsetbits_common(bitfield2));
    ASSERT_EQ(enumerator.counter(), bitfield2.nsetbits_common(bitfield1));
}
*/