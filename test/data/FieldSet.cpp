//
// Created by rja on 05/03/2020.
//
#if 0

#include <src/data/TableNew.h>
#include <src/data/FermionOnv.h>
#include "gtest/gtest.h"

TEST(FieldSet, Segments) {
    struct TestFieldSet : TableNew {
        Field<int> some_ints;

        TestFieldSet(size_t nsegment = 1) :
                TableNew(nsegment),
                some_ints(this, 17) {}
    };

    TestFieldSet table(4);
    table.extend(2);
}


TEST(FieldSet, Offsets) {
    struct TestTable : TableNew {
        Field<long int> long_ints;
        Field<short int> short_ints;
        Flag<> flaggg;
        DeterminantField<> det;

        TestTable() : TableNew(),
                      long_ints(this, 6),
                      short_ints(this, {7}),
                      flaggg(this, {6}),
                      det(this, 40) {}
    };
    TestTable table;

    table.extend(5);
    auto t = table.det(4);

    t[0].set(3);
    t[0].set(30);
    t[1].set(32);
    t[1].set(33);
    t.set(40);

    t.print();

    /*
    *field_set.short_ints(0, {0}) = 2;
    field_set.flaggg.set(0, {0});
    ASSERT_TRUE(field_set.flaggg.get(0, {0}));
    field_set.flaggg.clr(0, {0});
    ASSERT_FALSE(field_set.flaggg.get(0, {0}));

    ASSERT_EQ(field_set.m_ndatawords, 8);
    auto lengths = defs::inds{6, 7};
    ASSERT_EQ(field_set.field_nelements(), lengths);
    auto offsets = defs::inds{0, 6*sizeof(short int)};
    ASSERT_EQ(field_set.field_offsets(), lengths);
     */
}

TEST(FieldSet, NonNumeric){

    struct BitString: TableNew{
        BitStringField<> field;
        BitString(size_t nbit):TableNew(), field(this, nbit){}
    };

}
#endif