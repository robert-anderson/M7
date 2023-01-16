//
// Created by Robert J. Anderson on 09/07/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/sort/LocalExtremalRows.h"
#include <M7_lib/table/BufferedTable.h>

namespace local_extremal_rows_test {
    typedef SingleFieldRow<field::Number<int>> int_scalar_row_t;
    typedef buffered::Table<int_scalar_row_t> int_scalar_table_t;
    typedef SingleFieldRow<field::Number<std::complex<float>>> complex_scalar_row_t;
    typedef buffered::Table<complex_scalar_row_t> complex_scalar_table_t;
    static constexpr uint_t c_nfind = 4;

    static v_t<std::complex<float>> get_complex_data() {
        return {
                {1.0,  2.0},
                {1.5,  0.4},
                {1.1,  -2.1},
                {1.2,  -2.7},
                {5.1,  -2.9},
                {-3.2, 2.3},
                {-1.2, 2.5},
                {0.2,  -5.7}
        };
    }

    static void setup(int_scalar_table_t &table) {
        const v_t<int> data = {1, 4, -2, 6, 325, -234, 21, -7, 321};
        table.clear();
        auto row = table.m_row;
        for (auto &i: data) {
            row.push_back_jump();
            row.m_field = i;
        }
    }
    static void setup(complex_scalar_table_t &table) {
        auto data = get_complex_data();
        table.clear();
        auto row = table.m_row;
        for (auto &i: data) {
            row.push_back_jump();
            row.m_field = i;
        }
    }
}
TEST(LocalExtremalRows, EmptyTable) {
    using namespace local_extremal_rows_test;
    int_scalar_table_t table("test", {});
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<int> lxv(row1.m_field, row2.m_field, false, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), 0ul);
}

TEST(LocalExtremalRows, Ascending) {
    using namespace local_extremal_rows_test;
    int_scalar_table_t table("test", {});
    setup(table);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<int> lxv(row1.m_field, row2.m_field, false, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row1.jump(lxv[0]);
    ASSERT_EQ(row1.m_field, -234);
    row1.jump(lxv[1]);
    ASSERT_EQ(row1.m_field, -7);
    row1.jump(lxv[2]);
    ASSERT_EQ(row1.m_field, -2);
    row1.jump(lxv[3]);
    ASSERT_EQ(row1.m_field, 1);
}

TEST(LocalExtremalRows, AscendingAbs) {
    using namespace local_extremal_rows_test;
    int_scalar_table_t table("test", {});
    setup(table);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<int> lxv(row1.m_field, row2.m_field, false, true, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row1.jump(lxv[0]);
    ASSERT_EQ(row1.m_field, 1);
    row1.jump(lxv[1]);
    ASSERT_EQ(row1.m_field, -2);
    row1.jump(lxv[2]);
    ASSERT_EQ(row1.m_field, 4);
    row1.jump(lxv[3]);
    ASSERT_EQ(row1.m_field, 6);
}

TEST(LocalExtremalRows, Descending) {
    using namespace local_extremal_rows_test;
    int_scalar_table_t table("test", {});
    setup(table);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<int> lxv(row1.m_field, row2.m_field, true, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row1.jump(lxv[0]);
    ASSERT_EQ(row1.m_field, 325);
    row1.jump(lxv[1]);
    ASSERT_EQ(row1.m_field, 321);
    row1.jump(lxv[2]);
    ASSERT_EQ(row1.m_field, 21);
    row1.jump(lxv[3]);
    ASSERT_EQ(row1.m_field, 6);
}

TEST(LocalExtremalRows, DescendingAbs) {
    using namespace local_extremal_rows_test;
    int_scalar_table_t table("test", {});
    setup(table);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<int> lxv(row1.m_field, row2.m_field, true, true, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row1.jump(lxv[0]);
    ASSERT_EQ(row1.m_field, 325);
    row1.jump(lxv[1]);
    ASSERT_EQ(row1.m_field, 321);
    row1.jump(lxv[2]);
    ASSERT_EQ(row1.m_field, -234);
    row1.jump(lxv[3]);
    ASSERT_EQ(row1.m_field, 21);
}

TEST(LocalExtremalRows, AscendingComplex) {
    using namespace local_extremal_rows_test;
    complex_scalar_table_t table("test", {});
    setup(table);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<std::complex<float>> lxv(row1.m_field, row2.m_field, false, true, 0);
    lxv.find(c_nfind);
    // prepare verification data
    auto complex_data = get_complex_data();
    auto cmp_fn = [](std::complex<float> z1, std::complex<float> z2){return std::abs(z1) < std::abs(z2);};
    std::sort(complex_data.begin(), complex_data.end(), cmp_fn);
    // note that this is not a valid means of verification if any elements in the test data have the same magnitude
    for (uint_t irow = 0ul; irow < c_nfind; ++irow){
        row1.jump(lxv[irow]);
        const std::complex<float> &field = row1.m_field;
        ASSERT_EQ(field, complex_data[irow]);
    }
}

TEST(LocalExtremalRows, DescendingComplex) {
    using namespace local_extremal_rows_test;
    complex_scalar_table_t table("test", {});
    setup(table);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows<std::complex<float>> lxv(row1.m_field, row2.m_field, true, true, 0);
    lxv.find(c_nfind);
    // prepare verification data
    auto complex_data = get_complex_data();
    auto cmp_fn = [](std::complex<float> z1, std::complex<float> z2){return std::abs(z1) > std::abs(z2);};
    std::sort(complex_data.begin(), complex_data.end(), cmp_fn);
    // note that this is not a valid means of verification if any elements in the test data have the same magnitude
    for (uint_t irow = 0ul; irow < c_nfind; ++irow){
        row1.jump(lxv[irow]);
        const std::complex<float> &field = row1.m_field;
        ASSERT_EQ(field, complex_data[irow]);
    }
}

TEST(LocalExtremalRows, WithClearedRows) {
    using namespace local_extremal_rows_test;
    int_scalar_table_t table("test", {});
    setup(table);
    ASSERT_EQ(table.nrecord(), 9ul);
    ASSERT_EQ(table.nfreed_row(), 0ul);
    table.free(2); // -2
    ASSERT_EQ(table.nrecord(), 8ul);
    ASSERT_EQ(table.nfreed_row(), 1ul);
    table.free(7); // -7
    ASSERT_EQ(table.nrecord(), 7ul);
    ASSERT_EQ(table.nfreed_row(), 2ul);
    auto row1 = table.m_row;
    auto row2 = row1;
    LocalExtremalRows< int> lxv(row1.m_field, row2.m_field, false, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row1.jump(lxv[0]);
    ASSERT_EQ(row1.m_field, -234);
    row1.jump(lxv[1]);
    ASSERT_EQ(row1.m_field, 1);
    row1.jump(lxv[2]);
    ASSERT_EQ(row1.m_field, 4);
    row1.jump(lxv[3]);
    ASSERT_EQ(row1.m_field, 6);
}