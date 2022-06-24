//
// Created by Robert J. Anderson on 1/25/22.
//

#include "gtest/gtest.h"
#include "M7_lib/field/CompositeField.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/table/BufferedTable.h"


namespace composite_field_test {
    struct BitsAndInts : CompositeField<field::Bitset<uint16_t>, field::Numbers<int, 2>> {
        typedef CompositeField<field::Bitset<uint16_t>, field::Numbers<int, 2>> base_t;
        field::Bitset<uint16_t> m_bitset;
        field::Numbers<int, 2> m_ints;

        BitsAndInts(Row *row, uint_t nbit, uinta_t<2> ints_shape, std::string prefix) :
                base_t(m_bitset, m_ints),
                m_bitset(row, nbit, prefix + " bitset"),
                m_ints(row, ints_shape, prefix + " ints") {}

        BitsAndInts(const BitsAndInts &other) :
            base_t(m_bitset, m_ints), m_bitset(other.m_bitset), m_ints(other.m_ints) {}

        BitsAndInts &operator=(const BitsAndInts &other) {
            m_bitset = other.m_bitset;
            m_ints = other.m_ints;
            return *this;
        }
    };

    /**
     * composite fields are nestable to arbitrary depth
     */
    struct BitsAndIntsPair : CompositeField<BitsAndInts, BitsAndInts> {
        typedef CompositeField<BitsAndInts, BitsAndInts> base_t;
        BitsAndInts m_first;
        BitsAndInts m_second;
        BitsAndIntsPair(Row* row, uint_t nbit, uinta_t<2> ints_shape, std::string prefix):
            base_t(m_first, m_second),
            m_first(row, nbit, ints_shape, prefix + " first"),
            m_second(row, nbit, ints_shape, prefix + " second"){}

        BitsAndIntsPair(const BitsAndIntsPair& other):
            base_t(m_first, m_second), m_first(other.m_first), m_second(other.m_second){}

        BitsAndIntsPair& operator=(const BitsAndIntsPair& other) {
            m_first = other.m_first;
            m_second = other.m_second;
            return *this;
        }
    };

    struct TestRow : Row {
        BitsAndInts m_single;
        BitsAndIntsPair m_pair;
        TestRow(uint_t nbit, uinta_t<2> ints_shape):
            m_single(this, nbit, ints_shape, "single"), m_pair(this, nbit, ints_shape, "pair"){}
    };
}

TEST(CompositeField, DataIntegrity) {
    using namespace composite_field_test;
    const uint_t nbit = 9;
    const uinta_t<2> ints_shape = {4, 3};

    TestRow row(nbit, ints_shape);
    ASSERT_EQ(&row.m_single.m_bitset, &row.m_single.get<0>());
    ASSERT_EQ(row.m_single.m_bitset.m_format.m_nelement, nbit);
    ASSERT_EQ(&row.m_single.m_ints, &row.m_single.get<1>());
    ASSERT_EQ(row.m_single.m_ints.m_format.m_shape, ints_shape);
    ASSERT_EQ(&row.m_pair.m_first.m_bitset, &row.m_pair.get<0>().get<0>());
    ASSERT_EQ(&row.m_pair.m_first.m_ints, &row.m_pair.get<0>().get<1>());
    ASSERT_EQ(&row.m_pair.m_second.m_bitset, &row.m_pair.get<1>().get<0>());
    ASSERT_EQ(row.m_pair.m_second.m_bitset.m_format.m_nelement, nbit);
    ASSERT_EQ(&row.m_pair.m_second.m_ints, &row.m_pair.get<1>().get<1>());
    ASSERT_EQ(row.m_pair.m_second.m_ints.m_format.m_shape, ints_shape);

    ASSERT_EQ(row.nfield(), 6ul);
    for (uint_t i=0ul; i<row.nfield(); ++i) ASSERT_EQ(row.m_fields[i]->m_row, &row);
    ASSERT_EQ(row.m_child, nullptr);

    auto row_copy = row;
    for (uint_t i=0ul; i<row_copy.nfield(); ++i) ASSERT_EQ(row_copy.m_fields[i]->m_row, &row_copy);

    BufferedTable<TestRow> bt("test table", {{nbit, ints_shape}});
    bt.resize(10);
    bt.push_back(4);

    auto r1 = bt.m_row;
    r1.restart();
    r1.m_single.m_bitset = {0, 1, 3, 6, 7};
    r1.m_single.m_ints = {-1, 1, -2, 2, -3, 3, -4, -4, -5, 5, -6, 6};
    r1.m_pair.m_first.m_bitset = {0, 2, 4, 6};
    r1.m_pair.m_first.m_ints = {-6, 6, -5, 5, -4, 4, -3, -3, -2, 2, -1, 1};
    r1.m_pair.m_second.m_bitset = {5, 6, 8};
    r1.m_pair.m_second.m_ints = {-7, 7, -6, 6, -5, 5, -4, -4, -3, 3, -2, 2};
    auto r2 = bt.m_row;
    r2.jump(2);
    r2.m_single = r1.m_single;
    r2.m_pair = r1.m_pair;
    ASSERT_TRUE(r1.begin());
    ASSERT_TRUE(r2.begin());

    ASSERT_EQ(r1.m_single.m_bitset.m_row, &r1);
    ASSERT_EQ(r1.m_single.m_ints.m_row, &r1);
    ASSERT_EQ(r1.m_pair.m_first.m_ints.m_row, &r1);
    ASSERT_EQ(r1.m_pair.m_first.m_bitset.m_row, &r1);
    ASSERT_EQ(r1.m_pair.m_second.m_ints.m_row, &r1);
    ASSERT_EQ(r1.m_pair.m_second.m_bitset.m_row, &r1);

    ASSERT_EQ(r2.m_single.m_bitset.m_row, &r2);
    ASSERT_EQ(r2.m_single.m_ints.m_row, &r2);
    ASSERT_EQ(r2.m_pair.m_first.m_ints.m_row, &r2);
    ASSERT_EQ(r2.m_pair.m_first.m_bitset.m_row, &r2);
    ASSERT_EQ(r2.m_pair.m_second.m_ints.m_row, &r2);
    ASSERT_EQ(r2.m_pair.m_second.m_bitset.m_row, &r2);

    ASSERT_EQ(r2.m_single, r1.m_single);
    ASSERT_EQ(r2.m_pair, r1.m_pair);
}