//
// Created by rja on 02/10/2020.
//

#include "src/core/table/BufferedTable.h"
#include "src/core/field/Fields.h"
#include "gtest/gtest.h"


struct CommonFieldTypeTable : public TableX {
    fields::Numbers<short, 1> m_shorts;
    fields::Numbers<short, 1> m_more_shorts;

    CommonFieldTypeTable(size_t n1, size_t n2) :
            m_shorts(this, "some shorts", n1),
            m_more_shorts(this, "some more shorts", n2) {}
};

TEST(Table, CommonFieldOffset) {
    CommonFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.m_shorts.m_field.m_offset, 0);
    ASSERT_EQ(t1.m_shorts.m_field.m_size, 3 * sizeof(short));
    /*
     * check that the next field follows on gaplessly...
     */
    ASSERT_EQ(t1.m_more_shorts.m_field.m_offset, 3 * sizeof(short));
    ASSERT_EQ(t1.m_more_shorts.m_field.m_size, 4 * sizeof(short));
    /*
     * but not in this case, where the first field takes up a whole
     * number of datawords...
     */
    CommonFieldTypeTable t2(8, 5);
    ASSERT_EQ(t2.m_shorts.m_field.m_offset, 0);
    ASSERT_EQ(t2.m_shorts.m_field.m_size, 8 * sizeof(short));
    ASSERT_EQ(t2.m_more_shorts.m_field.m_offset, 2 * sizeof(defs::data_t));
    ASSERT_EQ(t2.m_more_shorts.m_field.m_size, 5 * sizeof(short));
}

struct DifferentFieldTypeTable : public TableX {
    fields::Numbers<short, 1> m_shorts;
    fields::Numbers<float, 1> m_floats;
    DifferentFieldTypeTable(size_t n1, size_t n2):
            m_shorts(this, "some shorts", n1),
            m_floats(this, "some floats", n2)
    {}
};

TEST(Table, DifferentFieldOffset){
    DifferentFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.m_shorts.m_field.m_offset, 0);
    ASSERT_EQ(t1.m_shorts.m_field.m_size, 3*sizeof(short));
    /*
     * check that the differently-typed second field is offset to the
     * next whole dataword
     */
    ASSERT_EQ(t1.m_floats.m_field.m_offset, sizeof(defs::data_t));
    ASSERT_EQ(t1.m_floats.m_field.m_size, 4*sizeof(float));

    ASSERT_EQ(t1.m_shorts.m_field.m_table, &t1);

    auto t2 = t1;
    ASSERT_EQ(t2.m_shorts.m_field.m_table, &t1);
}

struct TestFlagSet : FlagSet {
    Flag flag1;
    Flag flag2;
    Flags<1> flags1;
    Flags<1> flags2;
    TestFlagSet(fields::Bitset* bitset):
    FlagSet(bitset),
    flag1(this, "first flag"),
    flag2(this, "second flag"),
    flags1(this, "first rank-1 flag set", 6),
    flags2(this, "second rank-1 flag set", 6){}
};

struct FlagsTestTable : public TableX {
    fields::Flags<TestFlagSet> m_flags;
    FlagsTestTable(): m_flags(this, "Flagset"){}
};

TEST(Table, Flag){
    BufferedTable<FlagsTestTable> bt;
    bt.expand(1);
    bt.push_back();
    bt.print_field_details();
    bt.m_flags.flags1(0, 1) = 1;
    std::cout << bt.m_flags(0).to_string() << std::endl;
}

