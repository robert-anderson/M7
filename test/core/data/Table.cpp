//
// Created by rja on 02/10/2020.
//

#include "src/core/data/Table.h"
#include "src/core/data/Fields.h"
#include "src/core/data/NdCompositeFieldSingle.h"
//#include "src/core/data/.h"
#include "gtest/gtest.h"


struct CommonFieldTypeTable : public TableX {
    fields::Numbers<short, 1> m_shorts;
    fields::Numbers<short, 1> m_more_shorts;
    CommonFieldTypeTable(size_t n1, size_t n2):
            m_shorts(this, "some shorts", n1),
            m_more_shorts(this, "some more shorts", n2)
    {}
};

TEST(Table, CommonFieldOffset){
    CommonFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.m_shorts.offset(), 0);
    ASSERT_EQ(t1.m_shorts.size(), 3*sizeof(short));
    /*
     * check that the next field follows on gaplessly...
     */
    ASSERT_EQ(t1.m_more_shorts.offset(), 3*sizeof(short));
    ASSERT_EQ(t1.m_more_shorts.size(), 4*sizeof(short));
    /*
     * but not in this case, where the first field takes up a whole
     * number of datawords...
     */
    CommonFieldTypeTable t2(8, 5);
    ASSERT_EQ(t2.m_shorts.offset(), 0);
    ASSERT_EQ(t2.m_shorts.size(), 8*sizeof(short));
    ASSERT_EQ(t2.m_more_shorts.offset(), 2*sizeof(defs::data_t));
    ASSERT_EQ(t2.m_more_shorts.size(), 5*sizeof(short));
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
    ASSERT_EQ(t1.m_shorts.offset(), 0);
    ASSERT_EQ(t1.m_shorts.size(), 3*sizeof(short));
    /*
     * check that the differently-typed second field is offset to the
     * next whole dataword
     */
    ASSERT_EQ(t1.m_floats.offset(), sizeof(defs::data_t));
    ASSERT_EQ(t1.m_floats.size(), 4*sizeof(float));
}

namespace fieldsx {

    template <typename T, size_t nind>
    using Number = NdCompositeFieldSingle<NumericFieldX<T>, nind>;

}

struct Tt : TableX {
    fieldsx::Number<float, 1> m_num;
    Tt():m_num(this, {}, "ads", {4}){}
};

TEST(Table, TEST){
    Tt t;
}

//struct FlagsTestTable : public TableX {
//    FlagField<0> flag1;
//    FlagField<0> flag2;
//    FlagField<1> flags1;
//    FlagField<1> flags2;
//    FlagsTestTable():
//    flag1(this, "first flag"),
//    flag2(this, "second flag"),
//    flags1(this, "first rank-1 flag set", 6),
//    flags2(this, "second rank-1 flag set", 6)
//    {}
//};
//
//TEST(Table, Flag){
//    FlagsTestTable t;
//    t.print_field_details();
//}
