//
// Created by Robert J. Anderson on 02/10/2020.
//

#include <M7_lib/sample/PRNG.h>
#include <M7_lib/sort/ExtremalIndices.h>
#include <M7_lib/sort/LambdaQuickSorter.h>
#include "gtest/gtest.h"

#if 0
struct CommonFieldTypeTable : public Table {
    fields::Numbers<short, 1> m_shorts;
    fields::Numbers<short, 1> m_more_shorts;

    CommonFieldTypeTable(uint_t n1, uint_t n2) :
            m_shorts(this, "some shorts", n1),
            m_more_shorts(this, "some more shorts", n2) {}
};

TEST(Table, CommonFieldOffset) {
    CommonFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.m_shorts.m_column.m_offset, 0);
    ASSERT_EQ(t1.m_shorts.m_column.m_size, 3 * sizeof(short));
    /*
     * check that the next field follows on gaplessly...
     */
    ASSERT_EQ(t1.m_more_shorts.m_column.m_offset, 3 * sizeof(short));
    ASSERT_EQ(t1.m_more_shorts.m_column.m_size, 4 * sizeof(short));
    /*
     * but not in this case, where the first field takes up a whole
     * number of datawords...
     */
    CommonFieldTypeTable t2(8, 5);
    ASSERT_EQ(t2.m_shorts.m_column.m_offset, 0);
    ASSERT_EQ(t2.m_shorts.m_column.m_size, 8 * sizeof(short));
    ASSERT_EQ(t2.m_more_shorts.m_column.m_offset, 2 * sizeof(defs::data_t));
    ASSERT_EQ(t2.m_more_shorts.m_column.m_size, 5 * sizeof(short));
}

struct DifferentFieldTypeTable : public Table {
    fields::Numbers<short, 1> m_shorts;
    fields::Numbers<float, 1> m_floats;

    DifferentFieldTypeTable(uint_t n1, uint_t n2) :
            m_shorts(this, "some shorts", n1),
            m_floats(this, "some floats", n2) {}
};

TEST(Table, DifferentFieldOffset) {
    DifferentFieldTypeTable t1(3, 4);
    ASSERT_EQ(t1.m_shorts.m_column.m_offset, 0);
    ASSERT_EQ(t1.m_shorts.m_column.m_size, 3 * sizeof(short));
    /*
     * check that the differently-typed second field is offset to the
     * next whole dataword
     */
    ASSERT_EQ(t1.m_floats.m_column.m_offset, sizeof(defs::data_t));
    ASSERT_EQ(t1.m_floats.m_column.m_size, 4 * sizeof(float));

    ASSERT_EQ(t1.m_shorts.m_column.m_table, &t1);
}

struct TestFlagSet : FlagSet {
    Flag flag1;
    Flag flag2;
    Flags<1> flags1;
    Flags<1> flags2;

    TestFlagSet() :
            flag1(this, "first flag"),
            flag2(this, "second flag"),
            flags1(this, "first rank-1 flag set", 6),
            flags2(this, "second rank-1 flag set", 6) {}
};

struct FlagsTestTable : public Table {
    fields::Flags<TestFlagSet> m_flags;

    FlagsTestTable() : m_flags(this, "Flagset", {}) {}
};

TEST(Table, Flag) {
    BufferedTable<FlagsTestTable> bt("Flags test");
    bt.expand(1);
    bt.push_back();
    bt.print_column_details();
    bt.m_flags.flags1(0, 1) = 1;
    std::cout << bt.m_flags(0).to_string() << std::endl;
}


struct SortingTestTable : public Table {
    fields::Number<double> m_values;

    SortingTestTable() : m_values(this, "sortable values") {}
};

TEST(Table, Sorting) {
    BufferedTable<SortingTestTable> bt("Sorting test");
    const uint_t size = 100;
    PRNG prng(14, size);
    bt.expand(size);
    bt.push_back(size);
    for (uint_t i = 0ul; i < size; ++i) {
        bt.m_values(i) = prng.draw_float();
    }
    auto comp_fn = [&bt](const uint_t& i1, const uint_t i2){return bt.m_values(i1)>=bt.m_values(i2);};
    QuickSorter sorter(comp_fn);
    sorter.sort(bt);
    ASSERT_TRUE(sorter.is_sorted(bt));
}


TEST(Table, FieldBasedSorting) {
    BufferedTable<SortingTestTable> bt("Sorting test");
    const uint_t size = 100;
    PRNG prng(14, size);
    bt.expand(size);
    bt.push_back(size);
    for (uint_t i = 0ul; i < size; ++i) {
        bt.m_values(i) = prng.draw_float();
    }
    auto getter_fn = [&bt](const uint_t& i) -> const double& {return bt.m_values(i);};
    TableFieldSorter<fields::Number<double>> sorter(getter_fn);
    sorter.sort(bt);
    ASSERT_TRUE(sorter.is_sorted(bt));
}

TEST(Table, ExtremalIndices) {
    BufferedTable<SortingTestTable> bt("Extremal values test");
    const uint_t size = 100;
    PRNG prng(14, size);
    bt.expand(size);
    bt.push_back(size);
    for (uint_t i = 0ul; i < size; ++i) {
        bt.m_values(i) = prng.draw_float();
    }
    auto comp_fn = [&bt](const uint_t& i1, const uint_t i2){return bt.m_values(i1)<=bt.m_values(i2);};
    ExtremalIndices xv(comp_fn);
    xv.find(45);
    QuickSorter sorter(comp_fn);
    sorter.sort(bt);

    for (uint_t ifound=0ul; ifound<xv.nfound(); ++ifound){
        ASSERT_EQ(xv[ifound], sorter[ifound]);
    }
}


TEST(Table, FieldBasedExtremalValues) {
    BufferedTable<SortingTestTable> bt("Extremal values test");
    const uint_t size = 100;
    PRNG prng(14, size);
    bt.expand(size);
    bt.push_back(size);
    for (uint_t i = 0ul; i < size; ++i) {
        bt.m_values(i) = prng.draw_float();
    }

    TableExtremalValues<fields::Number<double>> xv(bt.m_values);
    xv.reset(bt);
    xv.find(45);
    TableFieldSorter<fields::Number<double>> sorter(bt.m_values);
    sorter.sort(bt);

    bt.print_contents(xv);

    for (uint_t ifound=0ul; ifound<xv.nfound(); ++ifound){
        ASSERT_EQ(xv[ifound], sorter[ifound]);
    }
}


TEST(Table, PointToPointTransfer) {
    if (mpi::nrank()==1) GTEST_SKIP();
    BufferedTable<DifferentFieldTypeTable> bt("Test", 3, 4);
    const uint_t ngen = 5;
    bt.push_back(ngen);
    for (uint_t i=0ul; i<ngen; ++i) bt.m_shorts(i, 0) = 55*(i+1)+mpi::irank();
    bt.print_contents();
    defs::uintv_t send;
    if (mpi::i_am(0)) send = {1, 2, 4};
    bt.transfer_rows(send, 0, 1);
    if (mpi::i_am(1)) {
        ASSERT_EQ(bt.m_hwm, 8);
    }
    std::cout << std::endl;
    bt.print_contents();
}

TEST(Table, Copy){
    const uint_t nshorts=4, nfloats=3;
    DifferentFieldTypeTable t(nshorts, nfloats);
    ASSERT_EQ(t.m_shorts.m_format.nelement(), nshorts);
    ASSERT_EQ(t.m_floats.m_format.nelement(), nfloats);
    auto tcopy = t;
    ASSERT_EQ(tcopy.m_last_copied, nullptr);
    ASSERT_EQ(t.m_last_copied, &tcopy);
    ASSERT_EQ(tcopy.m_row_size, t.m_row_size);
    ASSERT_EQ(tcopy.m_row_dsize, t.m_row_dsize);
    ASSERT_EQ(tcopy.m_shorts.m_format.nelement(), nshorts);
    ASSERT_EQ(tcopy.m_floats.m_format.nelement(), nfloats);
    ASSERT_EQ(tcopy.m_shorts.m_column.m_table, &tcopy);
    ASSERT_EQ(tcopy.m_floats.m_column.m_table, &tcopy);

    ASSERT_EQ(t.m_columns.size(), 2);
    ASSERT_EQ(tcopy.m_columns.size(), 2);
}

TEST(Table, CopyBuffered){
    const uint_t nshorts=4, nfloats=3;
    BufferedTable<DifferentFieldTypeTable> t("Copy test table", nshorts, nfloats);
    ASSERT_EQ(t.m_shorts.m_format.nelement(), nshorts);
    ASSERT_EQ(t.m_floats.m_format.nelement(), nfloats);

    t.push_back(4);
    t.m_shorts(0, 1) = 1023;
    t.m_floats(0, 0) = 10.23;
    t.m_shorts(1, 1) = 1123;
    t.m_floats(1, 0) = 11.23;
    t.m_shorts(2, 1) = 1223;
    t.m_floats(2, 0) = 12.23;
    t.m_shorts(3, 1) = 1323;
    t.m_floats(3, 0) = 13.23;

    auto tcopy = t;
    ASSERT_EQ(tcopy.m_last_copied, nullptr);
    ASSERT_EQ(t.m_last_copied, &tcopy);
    ASSERT_EQ(tcopy.m_row_size, t.m_row_size);
    ASSERT_EQ(tcopy.m_row_dsize, t.m_row_dsize);
    ASSERT_EQ(tcopy.m_shorts.m_format.nelement(), nshorts);
    ASSERT_EQ(tcopy.m_floats.m_format.nelement(), nfloats);
    ASSERT_EQ(tcopy.m_shorts.m_column.m_table, &tcopy);
    ASSERT_EQ(tcopy.m_floats.m_column.m_table, &tcopy);

    ASSERT_EQ(t.m_columns.size(), 2);
    ASSERT_EQ(tcopy.m_columns.size(), 2);

    ASSERT_EQ(t.m_hwm, 4);
    ASSERT_EQ(t.m_nrow, (uint_t)(4*(1+t.m_bw.expansion_factor())));

    // copied buffer should be empty
    ASSERT_EQ(tcopy.m_hwm, 0);
    ASSERT_EQ(tcopy.m_nrow, 0);
}
#endif