//
// Created by Robert J. Anderson on 07/06/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/foreach/Foreach.h"


/*
 * Compile-time number of dimensions
 */
/**
 * edge case where there are no dimensions to iterate over
 */
TEST(Foreach, CtndUnrestricted0) {
    const std::array<size_t, 0> shape = {};
    foreach::ctnd::Unrestricted<0> foreach(shape);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0);
    ASSERT_EQ(foreach.m_nterm, 0);
}

TEST(Foreach, CtndUnrestricted3) {
    const std::array<size_t, 3> shape = {3, 4, 2};
    const std::vector<std::array<size_t, 3>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {0, 2, 0},
            {0, 2, 1},
            {0, 3, 0},
            {0, 3, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1},
            {1, 2, 0},
            {1, 2, 1},
            {1, 3, 0},
            {1, 3, 1},
            {2, 0, 0},
            {2, 0, 1},
            {2, 1, 0},
            {2, 1, 1},
            {2, 2, 0},
            {2, 2, 1},
            {2, 3, 0},
            {2, 3, 1}
    };
    foreach::ctnd::Unrestricted<3> foreach(shape);
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, CtndOrderedStrictAsc0) {
    const size_t n = 5;
    foreach::ctnd::Ordered<0, true, true> foreach(n);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, CtndOrderedStrictAsc3) {
    const size_t n = 5;
    foreach::ctnd::Ordered<3, true, true> foreach(n);
    const std::vector<std::array<size_t, 3>> chk_inds = {
            {0, 1, 2},
            {0, 1, 3},
            {0, 2, 3},
            {1, 2, 3},
            {0, 1, 4},
            {0, 2, 4},
            {1, 2, 4},
            {0, 3, 4},
            {1, 3, 4},
            {2, 3, 4}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, CtndOrderedStrictDesc0) {
    const size_t n = 5;
    foreach::ctnd::Ordered<0, true, false> foreach(n);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, CtndOrderedStrictDesc3) {
    const size_t n = 5;
    foreach::ctnd::Ordered<3, true, false> foreach(n);
    const std::vector<std::array<size_t, 3>> chk_inds = {
            {2, 1, 0},
            {3, 1, 0},
            {3, 2, 0},
            {3, 2, 1},
            {4, 1, 0},
            {4, 2, 0},
            {4, 2, 1},
            {4, 3, 0},
            {4, 3, 1},
            {4, 3, 2}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, CtndOrderedAsc0) {
    const size_t n = 5;
    foreach::ctnd::Ordered<0, false, true> foreach(n);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, CtndOrderedAsc3) {
    const size_t n = 3;
    foreach::ctnd::Ordered<3, false, true> foreach(n);
    const std::vector<std::array<size_t, 3>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 1},
            {1, 1, 1},
            {0, 0, 2},
            {0, 1, 2},
            {1, 1, 2},
            {0, 2, 2},
            {1, 2, 2},
            {2, 2, 2}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, CtndOrderedDesc0) {
    const size_t n = 5;
    foreach::ctnd::Ordered<0, false, false> foreach(n);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, CtndOrderedDesc3) {
    const size_t n = 3;
    foreach::ctnd::Ordered<3, false, false> foreach(n);
    const std::vector<std::array<size_t, 3>> chk_inds = {
             {0, 0, 0},
             {1, 0, 0},
             {1, 1, 0},
             {1, 1, 1},
             {2, 0, 0},
             {2, 1, 0},
             {2, 1, 1},
             {2, 2, 0},
             {2, 2, 1},
             {2, 2, 2}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}


/*
 * Run-time number of dimensions
 */


/**
 * edge case where there are no dimensions to iterate over
 */
TEST(Foreach, RtndUnrestricted0) {
    const std::vector<size_t> shape = {};
    foreach::rtnd::Unrestricted foreach(shape);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0);
    ASSERT_EQ(foreach.m_nterm, 0);
}

TEST(Foreach, RtndUnrestricted3) {
    const std::vector<size_t> shape = {3, 4, 2};
    const std::vector<std::vector<size_t>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {0, 2, 0},
            {0, 2, 1},
            {0, 3, 0},
            {0, 3, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1},
            {1, 2, 0},
            {1, 2, 1},
            {1, 3, 0},
            {1, 3, 1},
            {2, 0, 0},
            {2, 0, 1},
            {2, 1, 0},
            {2, 1, 1},
            {2, 2, 0},
            {2, 2, 1},
            {2, 3, 0},
            {2, 3, 1}
    };
    foreach::rtnd::Unrestricted foreach(shape);
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, RtndOrderedStrictAsc0) {
    const size_t n = 5;
    foreach::rtnd::Ordered<true, true> foreach(n, 0);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, RtndOrderedStrictAsc3) {
    const size_t n = 5;
    foreach::rtnd::Ordered<true, true> foreach(n, 3);
    const std::vector<std::vector<size_t>> chk_inds = {
            {0, 1, 2},
            {0, 1, 3},
            {0, 2, 3},
            {1, 2, 3},
            {0, 1, 4},
            {0, 2, 4},
            {1, 2, 4},
            {0, 3, 4},
            {1, 3, 4},
            {2, 3, 4}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, RtndOrderedStrictDesc0) {
    const size_t n = 5;
    foreach::rtnd::Ordered<true, false> foreach(n, 0);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, RtndOrderedStrictDesc3) {
    const size_t n = 5;
    foreach::rtnd::Ordered<true, false> foreach(n, 3);
    const std::vector<std::vector<size_t>> chk_inds = {
            {2, 1, 0},
            {3, 1, 0},
            {3, 2, 0},
            {3, 2, 1},
            {4, 1, 0},
            {4, 2, 0},
            {4, 2, 1},
            {4, 3, 0},
            {4, 3, 1},
            {4, 3, 2}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, RtndOrderedAsc0) {
    const size_t n = 5;
    foreach::rtnd::Ordered<false, true> foreach(n, 0);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, RtndOrderedAsc3) {
    const size_t n = 3;
    foreach::rtnd::Ordered<false, true> foreach(n, 3);
    const std::vector<std::vector<size_t>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 1},
            {1, 1, 1},
            {0, 0, 2},
            {0, 1, 2},
            {1, 1, 2},
            {0, 2, 2},
            {1, 2, 2},
            {2, 2, 2}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}

TEST(Foreach, RtndOrderedDesc0) {
    const size_t n = 5;
    foreach::rtnd::Ordered<false, false> foreach(n, 0);
    size_t i = 0ul;
    auto fn = [&]() {++i;};
    foreach(fn);
    ASSERT_EQ(i, 0ul);
    ASSERT_EQ(foreach.m_nterm, 0ul);
}

TEST(Foreach, RtndOrderedDesc3) {
    const size_t n = 3;
    foreach::rtnd::Ordered<false, false> foreach(n, 3);
    const std::vector<std::vector<size_t>> chk_inds = {
            {0, 0, 0},
            {1, 0, 0},
            {1, 1, 0},
            {1, 1, 1},
            {2, 0, 0},
            {2, 1, 0},
            {2, 1, 1},
            {2, 2, 0},
            {2, 2, 1},
            {2, 2, 2}
    };
    size_t i = 0ul;
    auto fn = [&]() {ASSERT_EQ(foreach.inds(), chk_inds[i++]);};
    foreach(fn);
    ASSERT_EQ(i, foreach.m_nterm);
}