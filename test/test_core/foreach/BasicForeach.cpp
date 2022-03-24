//
// Created by rja on 24/03/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/foreach/BasicForeach.h"

/*
 * Compile-time number of dimensions
 */
/*
 * edge case where there are no dimensions to iterate over
 */
TEST(BasicForeach, CtndUnrestricted0) {
    using namespace basic_foreach;
    ctnd::Unrestricted<0> foreach({});
    ASSERT_EQ(foreach.m_niter, 0);
    size_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0> &inds) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndUnrestricted3) {
    using namespace basic_foreach;
    constexpr size_t nind = 3;
    const ctnd::inds_t<nind> shape = {3, 4, 2};
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
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
    ctnd::Unrestricted<nind> foreach(shape);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const ctnd::inds_t<nind>& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, CtndUnrestrictedExit) {
    using namespace basic_foreach;
    constexpr size_t nind = 3;
    const ctnd::inds_t<nind> shape = {2, 2, 2};
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1}
    };

    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        const auto& terminal_inds = chk_inds[i];

        size_t iiter = 0ul;
        auto fn = [&](const ctnd::inds_t<nind>& inds){
            if (inds==terminal_inds) throw ExitLoop();
            ++iiter;
        };

        ctnd::Unrestricted<nind> foreach(shape);
        foreach.loop(fn);
        // should have terminated on item i
        ASSERT_EQ(iiter, i);
    }
}

TEST(BasicForeach, CtndOrderedStrictAsc0) {
    using namespace basic_foreach;
    ctnd::Ordered<0, true, true> foreach(0);
    ASSERT_EQ(foreach.m_niter, 0);
    size_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0> &inds) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedStrictAsc3) {
    using namespace basic_foreach;
    constexpr size_t nind = 3;
    const size_t n = 5;
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
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
    ctnd::Ordered<nind, true, true> foreach(n);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const ctnd::inds_t<nind>& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, CtndOrderedExit) {
    using namespace basic_foreach;
    const size_t n = 5;
    constexpr size_t nind = 3;
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
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

    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        const auto& terminal_inds = chk_inds[i];

        size_t iiter = 0ul;
        auto fn = [&](const ctnd::inds_t<nind>& inds){
            if (inds==terminal_inds) throw ExitLoop();
            ++iiter;
        };

        ctnd::Ordered<nind, true, true> foreach(n);
        foreach.loop(fn);
        // should have terminated on item i
        ASSERT_EQ(iiter, i);
    }
}

TEST(BasicForeach, CtndOrderedStrictDesc0) {
    using namespace basic_foreach;
    ctnd::Ordered<0, true, false> foreach(0);
    ASSERT_EQ(foreach.m_niter, 0);
    size_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0> &inds) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedStrictDesc3) {
    using namespace basic_foreach;
    constexpr size_t nind = 3;
    const size_t n = 5;
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
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
    ctnd::Ordered<nind, true, false> foreach(n);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const ctnd::inds_t<nind>& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, CtndOrderedAsc0) {
    using namespace basic_foreach;
    ctnd::Ordered<0, false, true> foreach(0);
    ASSERT_EQ(foreach.m_niter, 0);
    size_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0> &inds) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedAsc3) {
    using namespace basic_foreach;
    constexpr size_t nind = 3;
    const size_t n = 3;
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
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
    ctnd::Ordered<nind, false, true> foreach(n);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const ctnd::inds_t<nind>& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, CtndOrderedDesc0) {
    using namespace basic_foreach;
    ctnd::Ordered<0, false, false> foreach(0);
    ASSERT_EQ(foreach.m_niter, 0);
    size_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0> &inds) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedDesc3) {
    using namespace basic_foreach;
    constexpr size_t nind = 3;
    const size_t n = 3;
    const std::vector<ctnd::inds_t<nind>> chk_inds = {
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
    ctnd::Ordered<nind, false, false> foreach(n);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const ctnd::inds_t<nind> &inds) {
        if (inds != *it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}