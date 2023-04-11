//
// Created by Robert J. Anderson on 24/03/2022.
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
    uint_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0>&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndUnrestricted3) {
    using namespace basic_foreach;
    constexpr uint_t nind = 3;
    const ctnd::inds_t<nind> shape = {3, 4, 2};
    const v_t<ctnd::inds_t<nind>> chk_inds = {
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
    constexpr uint_t nind = 3;
    const ctnd::inds_t<nind> shape = {2, 2, 2};
    const v_t<ctnd::inds_t<nind>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1}
    };

    for (uint_t i = 0ul; i < chk_inds.size(); ++i) {
        const auto& terminal_inds = chk_inds[i];

        uint_t iiter = 0ul;
        auto fn = [&](const ctnd::inds_t<nind>& inds){
            if (inds==terminal_inds) throw ExitLoop();
            ++iiter;
        };

        ctnd::Unrestricted<nind> foreach(shape);
        try { ASSERT_ANY_THROW(foreach.loop(fn)); }
        catch (const ExitLoop &) {}
        // should have terminated on item i
        ASSERT_EQ(iiter, i);
    }
}

TEST(BasicForeach, CtndOrderedStrictAsc0) {
    using namespace basic_foreach;
    ctnd::Ordered<0, true, true> foreach(0);
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0>&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedStrictAsc3) {
    using namespace basic_foreach;
    constexpr uint_t nind = 3;
    const uint_t n = 5;
    const v_t<ctnd::inds_t<nind>> chk_inds = {
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
    const uint_t n = 5;
    constexpr uint_t nind = 3;
    const v_t<ctnd::inds_t<nind>> chk_inds = {
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

    for (uint_t i = 0ul; i < chk_inds.size(); ++i) {
        const auto& terminal_inds = chk_inds[i];

        uint_t iiter = 0ul;
        auto fn = [&](const ctnd::inds_t<nind>& inds){
            if (inds==terminal_inds) throw ExitLoop();
            ++iiter;
        };

        ctnd::Ordered<nind, true, true> foreach(n);
        try {ASSERT_ANY_THROW(foreach.loop(fn));}
        catch (const ExitLoop &) {}

        // should have terminated on item i
        ASSERT_EQ(iiter, i);
    }
}

TEST(BasicForeach, CtndOrderedStrictDesc0) {
    using namespace basic_foreach;
    ctnd::Ordered<0, true, false> foreach(0);
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0>&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedStrictDesc3) {
    using namespace basic_foreach;
    constexpr uint_t nind = 3;
    const uint_t n = 5;
    const v_t<ctnd::inds_t<nind>> chk_inds = {
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
    uint_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0>&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedAsc3) {
    using namespace basic_foreach;
    constexpr uint_t nind = 3;
    const uint_t n = 3;
    const v_t<ctnd::inds_t<nind>> chk_inds = {
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
    uint_t iiter = 0ul;
    auto fn = [&](const ctnd::inds_t<0>&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, CtndOrderedDesc3) {
    using namespace basic_foreach;
    constexpr uint_t nind = 3;
    const uint_t n = 3;
    const v_t<ctnd::inds_t<nind>> chk_inds = {
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
    auto fn = [&](const ctnd::inds_t<nind>& inds) {
        if (inds != *it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

/*
 * Run-time number of dimensions
 */
TEST(BasicForeach, RtndUnrestricted0) {
    using namespace basic_foreach;
    rtnd::Unrestricted foreach({});
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const rtnd::inds_t&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, RtndUnrestricted3) {
    using namespace basic_foreach;
    const rtnd::inds_t shape = {3, 4, 2};
    const v_t<rtnd::inds_t> chk_inds = {
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
    rtnd::Unrestricted foreach(shape);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const rtnd::inds_t& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}


TEST(BasicForeach, RtndUnrestrictedExit) {
    using namespace basic_foreach;
    const rtnd::inds_t shape = {2, 2, 2};
    const v_t<rtnd::inds_t> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1}
    };

    for (uint_t i = 0ul; i < chk_inds.size(); ++i) {
        const auto& terminal_inds = chk_inds[i];

        uint_t iiter = 0ul;
        auto fn = [&](const rtnd::inds_t& inds){
            if (inds==terminal_inds) throw ExitLoop();
            ++iiter;
        };

        rtnd::Unrestricted foreach(shape);
        try { ASSERT_ANY_THROW(foreach.loop(fn)); }
        catch (const ExitLoop&){}
        // should have terminated on item i
        ASSERT_EQ(iiter, i);
    }
}


TEST(BasicForeach, RtndOrderedStrictAsc0) {
    using namespace basic_foreach;
    rtnd::Ordered<true, true> foreach(0, 0);
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const rtnd::inds_t&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, RtndOrderedStrictAsc3) {
    using namespace basic_foreach;
    const uint_t nind = 3;
    const uint_t n = 5;
    const v_t<rtnd::inds_t> chk_inds = {
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
    rtnd::Ordered<true, true> foreach(n, nind);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const rtnd::inds_t& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, RtndOrderedExit) {
    using namespace basic_foreach;
    const uint_t n = 5;
    const uint_t nind = 3;
    const v_t<rtnd::inds_t> chk_inds = {
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

    for (uint_t i = 0ul; i < chk_inds.size(); ++i) {
        const auto& terminal_inds = chk_inds[i];

        uint_t iiter = 0ul;
        auto fn = [&](const rtnd::inds_t& inds){
            if (inds==terminal_inds) throw ExitLoop();
            ++iiter;
        };

        rtnd::Ordered<true, true> foreach(n, nind);
        try { ASSERT_ANY_THROW(foreach.loop(fn)); }
        catch (const ExitLoop&){}
        // should have terminated on item i
        ASSERT_EQ(iiter, i);
    }
}

TEST(BasicForeach, RtndOrderedStrictDesc0) {
    using namespace basic_foreach;
    rtnd::Ordered<true, false> foreach(0, 0);
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const rtnd::inds_t&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, RtndOrderedStrictDesc3) {
    using namespace basic_foreach;
    const uint_t nind = 3;
    const uint_t n = 5;
    const v_t<rtnd::inds_t> chk_inds = {
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
    rtnd::Ordered<true, false> foreach(n, nind);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const rtnd::inds_t& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, RtndOrderedAsc0) {
    using namespace basic_foreach;
    rtnd::Ordered<false, true> foreach(0, 0);
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const rtnd::inds_t&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, RtndOrderedAsc3) {
    using namespace basic_foreach;
    const uint_t nind = 3;
    const uint_t n = 3;
    const v_t<rtnd::inds_t> chk_inds = {
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
    rtnd::Ordered<false, true> foreach(n, nind);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const rtnd::inds_t& inds) {
        if (inds!=*it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}

TEST(BasicForeach, RtndOrderedDesc0) {
    using namespace basic_foreach;
    rtnd::Ordered<false, false> foreach(0, 0);
    ASSERT_EQ(foreach.m_niter, 0);
    uint_t iiter = 0ul;
    auto fn = [&](const rtnd::inds_t&) { ++iiter; };
    foreach.loop(fn);
    ASSERT_FALSE(iiter);
}

TEST(BasicForeach, RtndOrderedDesc3) {
    using namespace basic_foreach;
    const uint_t nind = 3;
    const uint_t n = 3;
    const v_t<rtnd::inds_t> chk_inds = {
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
    rtnd::Ordered<false, false> foreach(n, nind);
    auto it = chk_inds.cbegin();
    bool all_correct = true;
    auto fn = [&](const rtnd::inds_t& inds) {
        if (inds != *it++) all_correct = false;
    };
    foreach.loop(fn);
    ASSERT_EQ(it, chk_inds.cend());
    ASSERT_TRUE(all_correct);
}