//
// Created by rja on 09/07/2021.
//

#include "gtest/gtest.h"
#include "src/core/sort/LocalExtremalValues.h"
#include <src/core/table/BufferedFields.h>

namespace local_extremal_values_test {
    typedef SingleFieldRow<fields::Number<int>> int_scalar_row_t;
    typedef BufferedTable<int_scalar_row_t> int_scalar_table_t;
    typedef SingleFieldRow<fields::Number<std::complex<float>>> complex_scalar_row_t;
    typedef BufferedTable<complex_scalar_row_t> complex_scalar_table_t;
    static constexpr size_t c_nfind = 4;

    static std::vector<std::complex<float>> get_complex_data() {
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
        const std::vector<int> data = {1, 4, -2, 6, 325, -234, 21, -7, 321};
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
TEST(LocalExtremalValues, EmptyTable) {
    using namespace local_extremal_values_test;
    int_scalar_table_t table("test", {{}});
    auto &row = table.m_row;
    LocalExtremalValues<int_scalar_row_t, int> lxv(row, row.m_field, false, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), 0ul);
}

TEST(LocalExtremalValues, Ascending) {
    using namespace local_extremal_values_test;
    int_scalar_table_t table("test", {{}});
    setup(table);
    auto &row = table.m_row;
    LocalExtremalValues<int_scalar_row_t, int> lxv(row, row.m_field, false, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row.jump(lxv[0]);
    ASSERT_EQ(row.m_field, -234);
    row.jump(lxv[1]);
    ASSERT_EQ(row.m_field, -7);
    row.jump(lxv[2]);
    ASSERT_EQ(row.m_field, -2);
    row.jump(lxv[3]);
    ASSERT_EQ(row.m_field, 1);
}

TEST(LocalExtremalValues, AscendingAbs) {
    using namespace local_extremal_values_test;
    int_scalar_table_t table("test", {{}});
    setup(table);
    auto &row = table.m_row;
    LocalExtremalValues<int_scalar_row_t, int> lxv(row, row.m_field, false, true, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row.jump(lxv[0]);
    ASSERT_EQ(row.m_field, 1);
    row.jump(lxv[1]);
    ASSERT_EQ(row.m_field, -2);
    row.jump(lxv[2]);
    ASSERT_EQ(row.m_field, 4);
    row.jump(lxv[3]);
    ASSERT_EQ(row.m_field, 6);
}

TEST(LocalExtremalValues, Descending) {
    using namespace local_extremal_values_test;
    int_scalar_table_t table("test", {{}});
    setup(table);
    auto &row = table.m_row;
    LocalExtremalValues<int_scalar_row_t, int> lxv(row, row.m_field, true, false, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row.jump(lxv[0]);
    ASSERT_EQ(row.m_field, 325);
    row.jump(lxv[1]);
    ASSERT_EQ(row.m_field, 321);
    row.jump(lxv[2]);
    ASSERT_EQ(row.m_field, 21);
    row.jump(lxv[3]);
    ASSERT_EQ(row.m_field, 6);
}

TEST(LocalExtremalValues, DescendingAbs) {
    using namespace local_extremal_values_test;
    int_scalar_table_t table("test", {{}});
    setup(table);
    auto &row = table.m_row;
    LocalExtremalValues<int_scalar_row_t, int> lxv(row, row.m_field, true, true, 0);
    lxv.find(c_nfind);
    ASSERT_EQ(lxv.nfound(), c_nfind);
    row.jump(lxv[0]);
    ASSERT_EQ(row.m_field, 325);
    row.jump(lxv[1]);
    ASSERT_EQ(row.m_field, 321);
    row.jump(lxv[2]);
    ASSERT_EQ(row.m_field, -234);
    row.jump(lxv[3]);
    ASSERT_EQ(row.m_field, 21);
}

TEST(LocalExtremalValues, AscendingComplex) {
    using namespace local_extremal_values_test;
    complex_scalar_table_t table("test", {{}});
    setup(table);
    auto &row = table.m_row;
    LocalExtremalValues<complex_scalar_row_t, std::complex<float>> lxv(row, row.m_field, false, true, 0);
    lxv.find(c_nfind);
    // prepare verification data
    auto complex_data = get_complex_data();
    auto cmp_fn = [](std::complex<float> z1, std::complex<float> z2){return std::abs(z1) < std::abs(z2);};
    std::sort(complex_data.begin(), complex_data.end(), cmp_fn);
    // note that this is not a valid means of verification if any elements in the test data have the same magnitude
    for (size_t irow = 0ul; irow < c_nfind; ++irow){
        row.jump(lxv[irow]);
        const std::complex<float> &field = row.m_field;
        ASSERT_EQ(field, complex_data[irow]);
    }
}

TEST(LocalExtremalValues, DescendingComplex) {
    using namespace local_extremal_values_test;
    complex_scalar_table_t table("test", {{}});
    setup(table);
    auto &row = table.m_row;
    LocalExtremalValues<complex_scalar_row_t, std::complex<float>> lxv(row, row.m_field, true, true, 0);
    lxv.find(c_nfind);
    // prepare verification data
    auto complex_data = get_complex_data();
    auto cmp_fn = [](std::complex<float> z1, std::complex<float> z2){return std::abs(z1) > std::abs(z2);};
    std::sort(complex_data.begin(), complex_data.end(), cmp_fn);
    // note that this is not a valid means of verification if any elements in the test data have the same magnitude
    for (size_t irow = 0ul; irow < c_nfind; ++irow){
        row.jump(lxv[irow]);
        const std::complex<float> &field = row.m_field;
        ASSERT_EQ(field, complex_data[irow]);
    }
}