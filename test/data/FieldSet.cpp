//
// Created by rja on 05/03/2020.
//

#include <src/data/FieldSet.h>
#include "gtest/gtest.h"

TEST(FieldSet, Offsets){
    struct TestFieldSet : FieldSet {
        Field<long int> long_ints;
        Field<short int> short_ints;
        TestFieldSet():FieldSet(nullptr),
        long_ints(this, 6),
        short_ints(this, 7)
        {}
    };
    TestFieldSet field_set;
    ASSERT_EQ(field_set.m_length, 8);
    auto lengths = defs::inds{6, 7};
    ASSERT_EQ(field_set.field_lengths(), lengths);
    auto offsets = defs::inds{0, 6*size_of(short int)};
    ASSERT_EQ(field_set.field_offsets(), lengths);
}
