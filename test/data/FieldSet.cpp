//
// Created by rja on 05/03/2020.
//

#include <src/data/FieldSet.h>
#include <src/data/Field.h>
#include "gtest/gtest.h"

TEST(FieldSet, Offsets){
    struct TestFieldSet : FieldSet {
        Field<long int> long_ints;
        Field<short int> short_ints;
        Flag<> flaggg;
        Bitfield<> bitfielddd;
        TestFieldSet(defs::inds &buffer):FieldSet(buffer.data()),
        long_ints(this, 6),
        short_ints(this, 7),
        flaggg(this, 6),
        bitfielddd(this, 65, 1)
        {}
    };
    defs::inds buffer(199);
    TestFieldSet field_set(buffer);
    *field_set.short_ints(0, 0) = 2;
    field_set.flaggg.set(0,0);
    ASSERT_TRUE(field_set.flaggg.get(0,0));
    field_set.flaggg.clr(0,0);
    ASSERT_FALSE(field_set.flaggg.get(0,0));

    /*
    ASSERT_EQ(field_set.m_length, 8);
    auto lengths = defs::inds{6, 7};
    ASSERT_EQ(field_set.field_nelements(), lengths);
    auto offsets = defs::inds{0, 6*sizeof(short int)};
    ASSERT_EQ(field_set.field_offsets(), lengths);
     */
}
