//
// Created by Robert J. Anderson on 07/06/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/foreach/ForeachVirtual.h"

#include <utility>

namespace foreach_virtual_test {
    using namespace foreach_virtual;

    template<size_t nind>
    struct CtndUnrestricted : ctnd::Unrestricted<nind> {
        using ctnd::Base<nind>::value;
        using ctnd::Base<nind>::iiter;
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<ctnd::inds_t<nind>> m_chk_inds;

        CtndUnrestricted(ctnd::inds_t<nind> shape, std::vector<ctnd::inds_t<nind>> chk_inds) :
                ctnd::Unrestricted<nind>(shape), m_chk_inds(std::move(chk_inds)) {}

        void body(const ctnd::inds_t<nind>& value, size_t iiter) override {
            DEBUG_ASSERT_LT(iiter, m_chk_inds.size(), "iteration count OOB");
            if (value != m_chk_inds[iiter]) m_pass = false;
        }
    };

    template<size_t nind, bool strict = true, bool ascending = true>
    struct CtndOrdered : ctnd::Ordered<nind, strict, ascending> {
        using ctnd::Base<nind>::value;
        using ctnd::Base<nind>::iiter;
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<ctnd::inds_t<nind>> m_chk_inds;

        CtndOrdered(size_t n, std::vector<ctnd::inds_t<nind>> chk_inds) :
                ctnd::Ordered<nind, strict, ascending>(n), m_chk_inds(std::move(chk_inds)) {}

        void body(const ctnd::inds_t<nind> &value, size_t iiter) override {
            DEBUG_ASSERT_LT(iiter, m_chk_inds.size(), "iteration count OOB");
            if (value != m_chk_inds[iiter]) m_pass = false;
        }
    };

    struct RtndUnrestricted : rtnd::Unrestricted {
        using rtnd::Base::value;
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<rtnd::inds_t> m_chk_inds;

        RtndUnrestricted(rtnd::inds_t shape, std::vector<rtnd::inds_t> chk_inds) :
                rtnd::Unrestricted(shape), m_chk_inds(std::move(chk_inds)) {}

        void body(const rtnd::inds_t &value, size_t iiter) override {
            DEBUG_ASSERT_LT(iiter, m_chk_inds.size(), "iteration count OOB");
            if (value!= m_chk_inds[iiter]) m_pass = false;
        }
    };

    template<bool strict = true, bool ascending = true>
    struct RtndOrdered : rtnd::Ordered<strict, ascending> {
        using rtnd::Base::value;
        using rtnd::Base::iiter;
        /**
         * pass status of the test
         */
        bool m_pass = true;
        /**
         * benchmark indices in order in which they should be generated by the loop
         */
        std::vector<rtnd::inds_t> m_chk_inds;

        RtndOrdered(size_t n, size_t r, std::vector<rtnd::inds_t> chk_inds) :
                rtnd::Ordered<strict, ascending>(n, r), m_chk_inds(std::move(chk_inds)) {}

        void body(const rtnd::inds_t &value, size_t iiter) override {
            DEBUG_ASSERT_LT(iiter, m_chk_inds.size(), "iteration count OOB");
            if (value != m_chk_inds[iiter]) m_pass = false;
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
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // testing lambda wrapper
    auto fn = [&chk_inds](const ctnd::inds_t<3>& value, size_t iiter){
        ASSERT_EQ(chk_inds[iiter], value);
    };
    ctnd::lambda::Unrestricted<3> lambda_foreach(fn, shape);
    lambda_foreach.loop();
    ASSERT_EQ(lambda_foreach.iiter()+1, lambda_foreach.niter());
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
        using CtndUnrestricted<3>::value;

        ctnd::inds_t<3> m_term_value;

        Foreach(ctnd::inds_t<3> shape, ctnd::inds_t<3> term_value) :
                CtndUnrestricted<3>(shape, {}), m_term_value(term_value) {}

        void body(const ctnd::inds_t<3> &value, size_t iiter) override {
            if (value== m_term_value) throw ExitLoop();
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(shape, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.iiter(), i);
    }
}


TEST(ForeachVirtual, CtndOrderedStrictAsc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, true, true> foreach(0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // testing lambda wrapper
    auto fn = [&chk_inds](const ctnd::inds_t<3>& value, size_t iiter){
        ASSERT_EQ(chk_inds[iiter], value);
    };
    ctnd::lambda::Ordered<3, true, true> lambda_foreach(fn, n);
    lambda_foreach.loop();
    ASSERT_EQ(lambda_foreach.iiter()+1, lambda_foreach.niter());
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
        ctnd::inds_t<3> m_term_value;

        Foreach(size_t n, ctnd::inds_t<3> term_value) :
                CtndOrdered<3, true, true>(n, {}), m_term_value(term_value) {}

        void body(const ctnd::inds_t<3> &value, size_t iiter) override {
            if (value == m_term_value) throw ExitLoop();
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(n, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.iiter(), i);
    }
}

TEST(ForeachVirtual, CtndOrderedStrictDesc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, true, false> foreach(0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedAsc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, false, true> foreach(0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, CtndOrderedDesc0) {
    using namespace foreach_virtual_test;
    CtndOrdered<0, false, false> foreach(0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);
}

/*
 * Run-time number of dimensions
 */
TEST(ForeachVirtual, RtndUnrestricted0) {
    using namespace foreach_virtual_test;
    RtndUnrestricted foreach({}, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // testing lambda wrapper
    auto fn = [&chk_inds](const defs::inds& value, size_t iiter){
        ASSERT_EQ(chk_inds[iiter], value);
    };
    rtnd::lambda::Unrestricted lambda_foreach(fn, shape);
    lambda_foreach.loop();
    ASSERT_EQ(lambda_foreach.iiter()+1, lambda_foreach.niter());
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
        rtnd::inds_t m_term_value;

        Foreach(rtnd::inds_t shape, rtnd::inds_t term_value) :
                RtndUnrestricted(std::move(shape), {}), m_term_value(std::move(term_value)) {}

        void body(const rtnd::inds_t &value, size_t iiter) override {
            if (value == m_term_value) throw ExitLoop();
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(shape, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.iiter(), i);
    }
}

TEST(ForeachVirtual, RtndOrderedStrictAsc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<true, true> foreach(0, 0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // testing lambda wrapper
    auto fn = [&chk_inds](const defs::inds& value, size_t iiter){
        ASSERT_EQ(chk_inds[iiter], value);
    };
    rtnd::lambda::Ordered<true, true> lambda_foreach(fn, n, 3);
    lambda_foreach.loop();
    ASSERT_EQ(lambda_foreach.iiter()+1, lambda_foreach.niter());
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
        rtnd::inds_t m_term_value;

        Foreach(size_t n, rtnd::inds_t term_value) :
                RtndOrdered<true, true>(n, 3, {}), m_term_value(std::move(term_value)) {}

        void body(const rtnd::inds_t &value, size_t iiter) override {
            if (value == m_term_value) throw ExitLoop();
        }
    };
    for (size_t i = 0ul; i < chk_inds.size(); ++i) {
        Foreach foreach(n, chk_inds[i]);
        foreach.loop();
        // should have terminated on item i
        ASSERT_EQ(foreach.iiter(), i);
    }
}

TEST(ForeachVirtual, RtndOrderedStrictDesc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<true, false> foreach(0, 0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedAsc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<false, true> foreach(0, 0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);
}

TEST(ForeachVirtual, RtndOrderedDesc0) {
    using namespace foreach_virtual_test;
    RtndOrdered<false, false> foreach(0, 0, {});
    ASSERT_EQ(foreach.niter(), 0);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 0);
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
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);

    // check the loop works as required when reused
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.niter());
    ASSERT_TRUE(foreach.m_pass);
}