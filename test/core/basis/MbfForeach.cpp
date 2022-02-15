//
// Created by anderson on 09/02/2022.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/basis/MbfForeach.h>
#include "gtest/gtest.h"

namespace mbf_foreach_test {
    namespace frm {
        namespace general {
            static std::vector<defs::inds> chk_inds() {
                return {
                        {0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}, {1, 2, 3, 4}, {0, 1, 2, 5},
                        {0, 1, 3, 5}, {0, 2, 3, 5}, {1, 2, 3, 5}, {0, 1, 4, 5}, {0, 2, 4, 5}, {1, 2, 4, 5},
                        {0, 3, 4, 5}, {1, 3, 4, 5}, {2, 3, 4, 5}
                };
            };
        }

        namespace spins {
            static std::vector<defs::inds> chk_inds() {
                return {
                        {0, 1, 6, 7}, {0, 2, 5, 7}, {1, 2, 4, 7}, {0, 3, 5, 6}, {1, 3, 4, 6}, {2, 3, 4, 5}
                };
            }
        }

        namespace ms2_conserve {
            static std::vector<defs::inds> chk_inds() {
                return {
                        {0, 1, 2, 4, 5}, {0, 1, 2, 4, 6}, {0, 1, 2, 5, 6}, {0, 1, 2, 4, 7}, {0, 1, 2, 5, 7},
                        {0, 1, 2, 6, 7},
                        {0, 1, 3, 4, 5}, {0, 1, 3, 4, 6}, {0, 1, 3, 5, 6}, {0, 1, 3, 4, 7}, {0, 1, 3, 5, 7},
                        {0, 1, 3, 6, 7},
                        {0, 2, 3, 4, 5}, {0, 2, 3, 4, 6}, {0, 2, 3, 5, 6}, {0, 2, 3, 4, 7}, {0, 2, 3, 5, 7},
                        {0, 2, 3, 6, 7},
                        {1, 2, 3, 4, 5}, {1, 2, 3, 4, 6}, {1, 2, 3, 5, 6}, {1, 2, 3, 4, 7}, {1, 2, 3, 5, 7},
                        {1, 2, 3, 6, 7}
                };
            }
        }
    }

    namespace bos {
        namespace general {
            static std::vector<defs::inds> chk_inds() {
                return {
                        {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 0}, {0, 1, 1}, {0, 1, 2}, {0, 2, 0}, {0, 2, 1},
                        {0, 2, 2},
                        {1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {1, 2, 0}, {1, 2, 1},
                        {1, 2, 2},
                        {2, 0, 0}, {2, 0, 1}, {2, 0, 2}, {2, 1, 0}, {2, 1, 1}, {2, 1, 2}, {2, 2, 0}, {2, 2, 1},
                        {2, 2, 2}
                };
            }
        }
    }
}


TEST(MbfForeach, FrmGeneral) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::general::chk_inds();
    auto fn = [&chk_inds](const field::FrmOnv &field, size_t iiter) {
        ASSERT_EQ(field, chk_inds[iiter]);
    };
    mbf_foreach::frm::General foreach(3, 4, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, chk_inds.size());
}

TEST(MbfForeach, FrmGeneralEarlyExit) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::general::chk_inds();
    auto fn = [&chk_inds](const field::FrmOnv &field, size_t iiter) {
        if (iiter==8) throw ExitLoop();
    };
    mbf_foreach::frm::General foreach(3, 4, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 8);
}

TEST(MbfForeach, FrmGeneralPair) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::general::chk_inds();
    auto fn = [&chk_inds](const field::FrmOnv &outer, size_t iouter, const field::FrmOnv &inner, size_t iinner) {
        ASSERT_EQ(outer, chk_inds[iouter]);
        ASSERT_EQ(inner, chk_inds[iinner]);
    };
    mbf_foreach::Pair<mbf_foreach::frm::General> foreach({3, 4, {}}, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, chk_inds.size() * chk_inds.size());
}

TEST(MbfForeach, FrmSpins) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::spins::chk_inds();
    auto fn = [&chk_inds](const field::FrmOnv &field, size_t iiter) {
        ASSERT_EQ(field, chk_inds[iiter]);
    };
    mbf_foreach::frm::Spins foreach(4, 0, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, chk_inds.size());
}

TEST(MbfForeach, FrmMs2Conserve) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::ms2_conserve::chk_inds();
    auto fn = [&chk_inds](const field::FrmOnv &field, size_t iiter) {
        ASSERT_EQ(field, chk_inds[iiter]);
    };
    mbf_foreach::frm::Ms2Conserve foreach(4, 5, 1, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, chk_inds.size());
}


TEST(MbfForeach, BosGeneral) {
    using namespace mbf_foreach_test;
    const auto chk_inds = bos::general::chk_inds();
    auto fn = [&chk_inds](const field::BosOnv &field, size_t iiter) {
        ASSERT_EQ(field, chk_inds[iiter]);
    };
    mbf_foreach::bos::General foreach(3, 2, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, chk_inds.size());
}

TEST(MbfForeach, BosGeneralEarlyExit) {
    using namespace mbf_foreach_test;
    auto fn = [](const field::BosOnv &field, size_t iiter) {
        if (iiter==8) throw ExitLoop();
    };
    mbf_foreach::bos::General foreach(3, 2, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 8);
}

TEST(MbfForeach, BosGeneralPair) {
    using namespace mbf_foreach_test;
    const auto chk_inds = bos::general::chk_inds();
    auto fn = [&chk_inds](const field::BosOnv &outer, size_t iouter, const field::BosOnv &inner, size_t iinner) {
        ASSERT_EQ(outer, chk_inds[iouter]);
        ASSERT_EQ(inner, chk_inds[iinner]);
    };
    mbf_foreach::Pair<mbf_foreach::bos::General> foreach({3, 2}, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, chk_inds.size() * chk_inds.size());
}

TEST(MbfForeach, FrmBosGeneral) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::general::chk_inds();
    const auto bos_chk_inds = bos::general::chk_inds();
    auto fn = [&frm_chk_inds, &bos_chk_inds](const field::FrmBosOnv &field, size_t iiter) {
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
    };
    mbf_foreach::frm::General outer(3, 4);
    mbf_foreach::bos::General inner(3, 2);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::General, mbf_foreach::bos::General> foreach(outer, inner, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, frm_chk_inds.size() * bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosSpins) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::spins::chk_inds();
    const auto bos_chk_inds = bos::general::chk_inds();
    auto fn = [&frm_chk_inds, &bos_chk_inds](const field::FrmBosOnv &field, size_t iiter) {
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
    };
    mbf_foreach::frm::Spins outer(4, 0);
    mbf_foreach::bos::General inner(3, 2);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::Spins, mbf_foreach::bos::General> foreach(outer, inner, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, frm_chk_inds.size() * bos_chk_inds.size());
}


TEST(MbfForeach, FrmBosSpinsEarlyExit) {
    using namespace mbf_foreach_test;
    auto fn = [](const field::FrmBosOnv &field, size_t iiter) {
        if (iiter==8) throw ExitLoop();
    };
    mbf_foreach::frm::Spins outer(4, 0);
    mbf_foreach::bos::General inner(3, 2);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::Spins, mbf_foreach::bos::General> foreach(outer, inner, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter(), 8);
}


TEST(MbfForeach, FrmBosMs2Conserve) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::ms2_conserve::chk_inds();
    const auto bos_chk_inds = bos::general::chk_inds();
    auto fn = [&frm_chk_inds, &bos_chk_inds](const field::FrmBosOnv &field, size_t iiter) {
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
    };
    mbf_foreach::frm::Ms2Conserve outer(4, 5, 1);
    mbf_foreach::bos::General inner(3, 2);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::Ms2Conserve, mbf_foreach::bos::General> foreach(outer, inner, fn);
    foreach.loop();
    ASSERT_EQ(foreach.iiter() + 1, frm_chk_inds.size() * bos_chk_inds.size());
}


TEST(MbfForeach, FrmBosMs2ConservePair) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::ms2_conserve::chk_inds();
    const auto bos_chk_inds = bos::general::chk_inds();
    auto fn = [&frm_chk_inds, &bos_chk_inds](const field::FrmBosOnv &outer, size_t iouter,
                                             const field::FrmBosOnv &inner, size_t iinner) {
        ASSERT_LT(iouter, frm_chk_inds.size() * bos_chk_inds.size());
        ASSERT_LT(iinner, frm_chk_inds.size() * bos_chk_inds.size());
        auto iouter_frm = iouter / bos_chk_inds.size();
        auto iouter_bos = iouter - iouter_frm * bos_chk_inds.size();
        auto iinner_frm = iinner / bos_chk_inds.size();
        auto iinner_bos = iinner - iinner_frm * bos_chk_inds.size();
        ASSERT_TRUE(outer.m_frm == frm_chk_inds[iouter_frm]);
        ASSERT_TRUE(outer.m_bos == bos_chk_inds[iouter_bos]);
        ASSERT_TRUE(inner.m_frm == frm_chk_inds[iinner_frm]);
        ASSERT_TRUE(inner.m_bos == bos_chk_inds[iinner_bos]);
    };
    mbf_foreach::frm::Ms2Conserve frm_foreach(4, 5, 1);
    mbf_foreach::bos::General bos_foreach(3, 2);
    typedef mbf_foreach::frm_bos::Product<mbf_foreach::frm::Ms2Conserve, mbf_foreach::bos::General> product_t;

    mbf_foreach::Pair<product_t> foreach(product_t(frm_foreach, bos_foreach), fn);
    foreach.loop();
    auto n = frm_chk_inds.size() * bos_chk_inds.size();
    n *= n;
    ASSERT_EQ(foreach.iiter() + 1, n);
}