//
// Created by Robert John Anderson on 2020-01-31.
//

#include <gtest/gtest.h>
/*
#include "../../src/data/Bitfield.h"
#include "../../src/enumerators/ClearedBitEnumerator.h"

TEST(ClearedBitEnumerator, Bitfield) {
    Bitfield bitfield1(200);
    Bitfield bitfield2(200);
    defs::inds tmp{};
    tmp.assign({1, 7, 8, 32, 89, 123, 199});
    bitfield1.set(tmp);
    tmp.assign({1, 4, 7, 32, 89, 125, 189, 190});
    bitfield2.set(tmp);
    ClearedBitEnumerator<true> enumerator;
    enumerator = ClearedBitEnumerator(bitfield1, bitfield2);

    tmp.assign({8, 123, 199});
    size_t setbit;
    while(enumerator.next(setbit)){
        ASSERT_TRUE(tmp[enumerator.counter()]==setbit);
    }
    ASSERT_EQ(enumerator.counter(), tmp.size());
    ASSERT_EQ(enumerator.counter(), bitfield1.nsetbits_cleared(bitfield2));

    // now the other way:
    enumerator = ClearedBitEnumerator(bitfield2, bitfield1);
    tmp.assign({4, 125, 189, 190});
    while(enumerator.next(setbit)){
        ASSERT_TRUE(tmp[enumerator.counter()]==setbit);
    }
    ASSERT_EQ(enumerator.counter(), tmp.size());
    ASSERT_EQ(enumerator.counter(), bitfield2.nsetbits_cleared(bitfield1));
}
*/