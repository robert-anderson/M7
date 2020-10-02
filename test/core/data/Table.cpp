//
// Created by rja on 02/10/2020.
//

#include "src/core/data/Table.h"
#include "src/core/data/NumericField.h"
#include "gtest/gtest.h"


struct CommonFieldTypeTable : public Table {
    NumericField<short, 1> shorts;
    NumericField<short, 1> more_shorts;
    CommonFieldTypeTable(size_t n1, size_t n2):
    shorts(this, n1),
    more_shorts(this, n2)
    {}
};

TEST(Table, CommonFieldOffset){
    CommonFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.shorts.offset(), 0);
    ASSERT_EQ(t1.shorts.size(), 3*sizeof(short));
    /*
     * check that the next field follows on gaplessly...
     */
    ASSERT_EQ(t1.more_shorts.offset(), 3*sizeof(short));
    ASSERT_EQ(t1.more_shorts.size(), 4*sizeof(short));
    /*
     * but not in this case, where the first field takes up a whole
     * number of datawords...
     */
    CommonFieldTypeTable t2(8, 5);
    ASSERT_EQ(t2.shorts.offset(), 0);
    ASSERT_EQ(t2.shorts.size(), 8*sizeof(short));
    ASSERT_EQ(t2.more_shorts.offset(), 2*sizeof(defs::data_t));
    ASSERT_EQ(t2.more_shorts.size(), 5*sizeof(short));
}

struct DifferentFieldTypeTable : public Table {
    NumericField<short, 1> shorts;
    NumericField<float, 1> floats;
    DifferentFieldTypeTable(size_t n1, size_t n2):
            shorts(this, n1),
            floats(this, n2)
    {}
};

TEST(Table, DifferentFieldOffset){
    DifferentFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.shorts.offset(), 0);
    ASSERT_EQ(t1.shorts.size(), 3*sizeof(short));
    /*
     * check that the differently-typed second field is offset to the
     * next whole dataword
     */
    ASSERT_EQ(t1.floats.offset(), sizeof(defs::data_t));
    ASSERT_EQ(t1.floats.size(), 4*sizeof(float));
}
