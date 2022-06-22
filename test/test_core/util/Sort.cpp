//
// Created by rja on 15/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Sort.h"

TEST(UtilSort, RealAsc) {
    const bool asc = true, abs_val = false;
    std::vector<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    defs::inds_t inds_chk = {5, 3, 0, 1, 2, 4};
    auto inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    std::vector<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    utils::sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // inds_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, RealDesc) {
    const bool asc = false, abs_val = false;
    std::vector<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    defs::inds_t inds_chk = {4, 2, 1, 0, 3, 5};
    auto inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    std::vector<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    utils::sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // inds_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, RealAscAbs) {
    const bool asc = true, abs_val = true;
    std::vector<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    defs::inds_t inds_chk = {0, 1, 3, 2, 4, 5};
    auto inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    std::vector<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    utils::sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // inds_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}

TEST(UtilSort, RealDescAbs) {
    const bool asc = false, abs_val = true;
    std::vector<double> v = {1.0, 2.3, 4.0, -2.3, 9.0, -10.3};
    defs::inds_t inds_chk = {5, 4, 2, 1, 3, 0};
    auto inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
    std::vector<double> v_chk;
    for (auto i: inds) v_chk.push_back(v[i]);
    utils::sort::inplace(v, asc, abs_val);
    ASSERT_EQ(v, v_chk);
    // inds_t should already be in order when called again after inplace sort
    std::iota(inds_chk.begin(), inds_chk.end(), 0ul);
    inds = utils::sort::inds(v, asc, abs_val);
    ASSERT_EQ(inds, inds_chk);
}