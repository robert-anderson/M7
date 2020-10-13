//
// Created by rja on 02/10/2020.
//

#include "src/core/data/Table_NEW.h"
#include "src/core/data/NumericField.h"
#include "src/core/data/FlagField.h"
#include "gtest/gtest.h"


struct CommonFieldTypeTable : public Table_NEW {
    NumericField<short, 1> shorts;
    NumericField<short, 1> more_shorts;
    CommonFieldTypeTable(size_t n1, size_t n2):
    shorts(this, "some shorts", n1),
    more_shorts(this, "some more shorts", n2)
    {}
};

TEST(Table, CommonFieldOffset){
    CommonFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.shorts.m_offset, 0);
    ASSERT_EQ(t1.shorts.m_size, 3*sizeof(short));
    /*
     * check that the next field follows on gaplessly...
     */
    ASSERT_EQ(t1.more_shorts.m_offset, 3*sizeof(short));
    ASSERT_EQ(t1.more_shorts.m_size, 4*sizeof(short));
    /*
     * but not in this case, where the first field takes up a whole
     * number of datawords...
     */
    CommonFieldTypeTable t2(8, 5);
    ASSERT_EQ(t2.shorts.m_offset, 0);
    ASSERT_EQ(t2.shorts.m_size, 8*sizeof(short));
    ASSERT_EQ(t2.more_shorts.m_offset, 2*sizeof(defs::data_t));
    ASSERT_EQ(t2.more_shorts.m_size, 5*sizeof(short));
}

struct DifferentFieldTypeTable : public Table_NEW {
    NumericField<short, 1> shorts;
    NumericField<float, 1> floats;
    DifferentFieldTypeTable(size_t n1, size_t n2):
            shorts(this, "some shorts", n1),
            floats(this, "some floats", n2)
    {}
};

TEST(Table, DifferentFieldOffset){
    DifferentFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.shorts.m_offset, 0);
    ASSERT_EQ(t1.shorts.m_size, 3*sizeof(short));
    /*
     * check that the differently-typed second field is offset to the
     * next whole dataword
     */
    ASSERT_EQ(t1.floats.m_offset, sizeof(defs::data_t));
    ASSERT_EQ(t1.floats.m_size, 4*sizeof(float));
}

struct FlagsTestTable : public Table_NEW {
    FlagField<0> flag1;
    FlagField<0> flag2;
    FlagField<1> flags1;
    FlagField<1> flags2;
    FlagsTestTable():
    flag1(this, "first flag"),
    flag2(this, "second flag"),
    flags1(this, "first rank-1 flag set", 6),
    flags2(this, "second rank-1 flag set", 6)
    {}
};

TEST(Table, Flag){
    FlagsTestTable t;
    t.print_field_details();
}
