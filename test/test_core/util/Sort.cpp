//
// Created by rja on 15/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Sort.h"

TEST(UtilSort, RealAsc) {
    const bool asc = true, abs_val = false;
    v_t<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    uintv_t inds_chk = {5, 3, 0, 1, 2, 4};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, RealDesc) {
    const bool asc = false, abs_val = false;
    v_t<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    uintv_t inds_chk = {4, 2, 1, 0, 3, 5};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, RealAscAbs) {
    const bool asc = true, abs_val = true;
    v_t<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    uintv_t inds_chk = {0, 1, 3, 2, 4, 5};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, RealDescAbs) {
    const bool asc = false, abs_val = true;
    v_t<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    uintv_t inds_chk = {5, 4, 2, 1, 3, 0};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, ComplexAsc) {
    const bool asc = true, abs_val = false;
    v_t<std::complex<double>> v = {{1.0, -0.1}, {2.3, -0.1}, {4.0, 0.1}, {-2.3, 0.1}, {9.0, -0.1}, {-10.3, 0.1}};
    uintv_t inds_chk = {5, 3, 0, 1, 2, 4};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<std::complex<double>> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, ComplexDesc) {
    const bool asc = false, abs_val = false;
    v_t<std::complex<double>> v = {{1.0, -0.1}, {2.3, -0.1}, {4.0, 0.1}, {-2.3, 0.1}, {9.0, -0.1}, {-10.3, 0.1}};
    uintv_t inds_chk = {4, 2, 1, 0, 3, 5};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<std::complex<double>> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, ComplexAscAbs) {
    const bool asc = true, abs_val = true;
    v_t<std::complex<double>> v = {{1.0, -0.1}, {2.3, -0.1}, {4.0, 0.1}, {-2.3, 0.1}, {9.0, -0.1}, {-10.3, 0.1}};
    uintv_t inds_chk = {0, 1, 3, 2, 4, 5};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<std::complex<double>> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, ComplexDescAbs) {
    const bool asc = false, abs_val = true;
    v_t<std::complex<double>> v = {{1.0, -0.1}, {2.3, -0.1}, {4.0, 0.1}, {-2.3, 0.1}, {9.0, -0.1}, {-10.3, 0.1}};
    uintv_t inds_chk = {5, 4, 2, 1, 3, 0};
    auto inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    v_t<std::complex<double>> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // uintv_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}
