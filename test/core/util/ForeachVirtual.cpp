//
// Created by rja on 07/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/util/ForeachVirtual.h"

#include <utility>

namespace foreach_virtual_test {
    using namespace foreach_virtual;

    template<size_t nind>
    struct CtndUnrestricted : ctnd::Unrestricted<nind> {
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<ctnd::inds_t<nind>> m_chk_inds;
        /**
         * body call counter
         */
        size_t m_i = 0ul;

        CtndUnrestricted(ctnd::inds_t<nind> shape, std::vector<ctnd::inds_t<nind>> chk_inds) :
                ctnd::Unrestricted<nind>(shape), m_chk_inds(std::move(chk_inds)) {}

        void body(const ctnd::inds_t<nind> &inds) override {
            if (inds != m_chk_inds[m_i++]) m_pass = false;
        }
    };

    template<size_t nind, bool strict = true, bool ascending = true>
    struct CtndOrdered : ctnd::Ordered<nind, strict, ascending> {
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<ctnd::inds_t<nind>> m_chk_inds;
        /**
         * body call counter
         */
        size_t m_i = 0ul;

        CtndOrdered(size_t n, std::vector<ctnd::inds_t<nind>> chk_inds) :
                ctnd::Ordered<nind, strict, ascending>(n), m_chk_inds(std::move(chk_inds)) {}

        void body(const ctnd::inds_t<nind> &inds) override {
            if (inds != m_chk_inds[m_i++]) m_pass = false;
        }
    };

    struct RtndUnrestricted : rtnd::Unrestricted {
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<rtnd::inds_t> m_chk_inds;
        /**
         * body call counter
         */
        size_t m_i = 0ul;

        RtndUnrestricted(rtnd::inds_t shape, std::vector<rtnd::inds_t> chk_inds) :
                rtnd::Unrestricted(shape), m_chk_inds(std::move(chk_inds)) {}

        void body(const rtnd::inds_t &inds) override {
            if (inds != m_chk_inds[m_i++]) m_pass = false;
        }
    };

    template<bool strict = true, bool ascending = true>
    struct RtndOrdered : rtnd::Ordered<strict, ascending> {
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<rtnd::inds_t> m_chk_inds;
        /**
         * body call counter
         */
        size_t m_i = 0ul;

        RtndOrdered(size_t n, size_t r, std::vector<rtnd::inds_t> chk_inds) :
                rtnd::Ordered<strict, ascending>(n, r), m_chk_inds(std::move(chk_inds)) {}

        void body(const rtnd::inds_t &inds) override {
            if (inds != m_chk_inds[m_i++]) m_pass = false;
        }
    };
}


/*
 * Compile-time number of dimensions
 */
/*
 * edge case where there are no dimensions to iterate over
 */
TEST(ForeachVirtual, CtndUnrestricted0) {
    using namespace foreach_virtual_test;
    CtndUnrestricted<0> foreach({}, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndUnrestricted3) {
    using namespace foreach_virtual_test;
    const ctnd::inds_t<3> shape = {3, 4, 2};
    const std::vector<ctnd::inds_t<3>> chk_inds = {
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
    CtndUnrestricted<3> foreach(shape, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndUnrestrictedExit) {
    using namespace foreach_virtual_test;
    const ctnd::inds_t<3> shape = {2, 2, 2};
    const std::vector<ctnd::inds_t<3>> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1},
    };
    struct Foreach : CtndUnrestricted<3> {
        ctnd::inds_t<3> m_term_inds;

        Foreach(ctnd::inds_t<3> shape, ctnd::inds_t<3> term_inds) :
                CtndUnrestricted<3>(shape, {}), m_term_inds(term_inds) {}

        void body(const ctnd::inds_t<3> &inds) override {
            if (inds == m_term_inds) throw ExitLoop();
            ++m_i;
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(shape, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.m_i, i);
    }
}

TEST(ForeachVirtual, CtndOrderedStrictAsc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, true, true> foreach(0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedStrictAsc3) {
    using namespace foreach_virtual_test;
    const size_t n = 5;
    const std::vector<ctnd::inds_t<3>> chk_inds = {
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
    CtndOrdered<3, true, true> foreach(n, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedExit) {
    using namespace foreach_virtual_test;
    const size_t n = 5;
    const std::vector<ctnd::inds_t<3>> chk_inds = {
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

    struct Foreach : CtndOrdered<3, true, true> {
        ctnd::inds_t<3> m_term_inds;

        Foreach(size_t n, ctnd::inds_t<3> term_inds) :
                CtndOrdered<3, true, true>(n, {}), m_term_inds(term_inds) {}

        void body(const ctnd::inds_t<3> &inds) override {
            if (inds == m_term_inds) throw ExitLoop();
            ++m_i;
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(n, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.m_i, i);
    }
}

TEST(ForeachVirtual, CtndOrderedStrictDesc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, true, false> foreach(0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedStrictDesc3) {
    using namespace foreach_virtual_test;
    const size_t n = 5;
    const std::vector<ctnd::inds_t<3>> chk_inds = {
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
    CtndOrdered<3, true, false> foreach(n, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedAsc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, false, true> foreach(0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedAsc3) {
    using namespace foreach_virtual_test;
    const size_t n = 3;
    const std::vector<ctnd::inds_t<3>> chk_inds = {
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
    CtndOrdered<3, false, true> foreach(n, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedDesc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, false, false> foreach(0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedDesc3) {
    using namespace foreach_virtual_test;
    const size_t n = 3;
    const std::vector<ctnd::inds_t<3>> chk_inds = {
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
    CtndOrdered<3, false, false> foreach(n, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

/*
 * Run-time number of dimensions
 */
TEST(ForeachVirtual, RtndUnrestricted0) {
    using namespace foreach_virtual_test;
    RtndUnrestricted foreach({}, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndUnrestricted3) {
    using namespace foreach_virtual_test;
    const rtnd::inds_t shape = {3, 4, 2};
    const std::vector<rtnd::inds_t> chk_inds = {
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
    RtndUnrestricted foreach(shape, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndUnrestrictedExit) {
    using namespace foreach_virtual_test;
    const rtnd::inds_t shape = {2, 2, 2};
    const std::vector<rtnd::inds_t> chk_inds = {
            {0, 0, 0},
            {0, 0, 1},
            {0, 1, 0},
            {0, 1, 1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, 0},
            {1, 1, 1},
    };
    struct Foreach : RtndUnrestricted {
        rtnd::inds_t m_term_inds;

        Foreach(rtnd::inds_t shape, rtnd::inds_t term_inds) :
                RtndUnrestricted(std::move(shape), {}), m_term_inds(std::move(term_inds)) {}

        void body(const rtnd::inds_t &inds) override {
            if (inds == m_term_inds) throw ExitLoop();
            ++m_i;
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(shape, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.m_i, i);
    }
}

TEST(ForeachVirtual, RtndOrderedStrictAsc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<true, true> foreach(0, 0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedStrictAsc3) {
    using namespace foreach_virtual_test;
    const size_t n = 5;
    const std::vector<rtnd::inds_t> chk_inds = {
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
    RtndOrdered<true, true> foreach(n, 3, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedExit) {
    using namespace foreach_virtual_test;
    const size_t n = 5;
    const std::vector<rtnd::inds_t> chk_inds = {
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

    struct Foreach : RtndOrdered<true, true> {
        rtnd::inds_t m_term_inds;

        Foreach(size_t n, rtnd::inds_t term_inds) :
                RtndOrdered<true, true>(n, 3, {}), m_term_inds(std::move(term_inds)) {}

        void body(const rtnd::inds_t &inds) override {
            if (inds == m_term_inds) throw ExitLoop();
            ++m_i;
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(n, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.m_i, i);
    }
}

TEST(ForeachVirtual, RtndOrderedStrictDesc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<true, false> foreach(0, 0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedStrictDesc3) {
    using namespace foreach_virtual_test;
    const size_t n = 5;
    const std::vector<rtnd::inds_t> chk_inds = {
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
    RtndOrdered<true, false> foreach(n, 3, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedAsc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<false, true> foreach(0, 0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedAsc3) {
    using namespace foreach_virtual_test;
    const size_t n = 3;
    const std::vector<rtnd::inds_t> chk_inds = {
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
    RtndOrdered<false, true> foreach(n, 3,chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedDesc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<false, false> foreach(0, 0, {});
    ASSERT_EQ(foreach.m_nterm, 0);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, 0);
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedDesc3) {
    using namespace foreach_virtual_test;
    const size_t n = 3;
    const std::vector<rtnd::inds_t> chk_inds = {
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
    RtndOrdered<false, false> foreach(n, 3, chk_inds);
    foreach.loop();
    ASSERT_EQ(foreach.m_i, foreach.m_nterm);
    ASSERT_TRUE(foreach.m_pass);
}