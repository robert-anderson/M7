//
// Created by Robert J. Anderson on 09/02/2022.
//

#include "gtest/gtest.h"

#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/foreach/MbfForeach.h"
#include "BasisData.h"

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
        namespace general_open {
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
    const sys::frm::Basis hs(4, 3);
    buffered::FrmOnv mbf(hs);
    size_t iiter = 0ul;
    auto fn = [&](const field::FrmOnv &mbf) {
        ASSERT_EQ(mbf, chk_inds[iiter]);
        ++iiter;
    };
    mbf_foreach::frm::General foreach(hs);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, chk_inds.size());
    ASSERT_EQ(foreach.m_niter, chk_inds.size());
}

TEST(MbfForeach, FrmGeneralEarlyExit) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::general::chk_inds();
    const sys::frm::Basis hs(4, 3);
    size_t iiter = 0ul;
    auto fn = [&](const field::FrmOnv &mbf) {
        ASSERT_EQ(mbf, chk_inds[iiter]);
        if (iiter == 8) throw ExitLoop();
        ++iiter;
    };
    mbf_foreach::frm::General foreach(hs);
    try { ASSERT_ANY_THROW(foreach.loop_fn(fn)); }
    catch (const ExitLoop &) {}
    ASSERT_EQ(iiter, 8);
    ASSERT_EQ(foreach.m_niter, chk_inds.size());
}

TEST(MbfForeach, FrmGeneralPair) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::general::chk_inds();
    const sys::frm::Basis hs(4, 3);
    buffered::FrmOnv mbf(hs);

    auto fn = [&chk_inds](
            const field::FrmOnv &outer, size_t iouter, const field::FrmOnv &inner, size_t iinner) {
        ASSERT_EQ(outer, chk_inds[iouter]);
        ASSERT_EQ(inner, chk_inds[iinner]);
    };
    mbf_foreach::frm::Pair<mbf_foreach::frm::General> foreach(hs);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(foreach.m_niter, chk_inds.size()*chk_inds.size());
}

TEST(MbfForeach, FrmSpins) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::spins::chk_inds();
    const sys::frm::Basis hs(4, 4, 0);
    buffered::FrmOnv mbf(hs);
    size_t iiter = 0ul;
    auto fn = [&](const field::FrmOnv &mbf) {
        ASSERT_EQ(mbf, chk_inds[iiter]);
        ++iiter;
    };
    mbf_foreach::frm::Spins foreach(hs);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, chk_inds.size());
    ASSERT_EQ(foreach.m_niter, chk_inds.size());
}

TEST(MbfForeach, FrmMs2Conserve) {
    using namespace mbf_foreach_test;
    const auto chk_inds = frm::ms2_conserve::chk_inds();
    const sys::frm::Basis hs(5, 4, 1);
    buffered::FrmOnv mbf(hs);
    size_t iiter = 0ul;
    auto fn = [&](const field::FrmOnv &mbf) {
        ASSERT_EQ(mbf, chk_inds[iiter]);
        ++iiter;
    };
    mbf_foreach::frm::Ms2Conserve foreach(hs);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, chk_inds.size());
    ASSERT_EQ(foreach.m_niter, chk_inds.size());
}

TEST(MbfForeach, BosGeneralOpen) {
    using namespace mbf_foreach_test;
    const auto chk_inds = bos::general_open::chk_inds();
    const BosHilbertSpace hs(0, 3, false, 2);
    buffered::BosOnv mbf(hs);
    size_t iiter = 0ul;
    auto fn = [&](const field::BosOnv &mbf) {
        ASSERT_EQ(mbf, chk_inds[iiter]);
        ++iiter;
    };
    mbf_foreach::bos::GeneralOpen foreach(hs);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, chk_inds.size());
    ASSERT_EQ(foreach.m_niter, chk_inds.size());
}

TEST(MbfForeach, BosGeneralOpenEarlyExit) {
    using namespace mbf_foreach_test;
    const auto chk_inds = bos::general_open::chk_inds();
    const BosHilbertSpace hs(0, 3, false, 2);
    size_t iiter = 0ul;
    auto fn = [&](const field::BosOnv &field) {
        ASSERT_EQ(field, chk_inds[iiter]);
        if (iiter == 8) throw ExitLoop();
        ++iiter;
    };
    mbf_foreach::bos::GeneralOpen foreach(hs);
    ASSERT_ANY_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, 8);
    ASSERT_EQ(foreach.m_niter, chk_inds.size());
}

TEST(MbfForeach, BosGeneralOpenPair) {
    using namespace mbf_foreach_test;
    const auto chk_inds = bos::general_open::chk_inds();
    const BosHilbertSpace hs(0, 3, false, 2);
    auto fn = [&chk_inds](
            const field::BosOnv &outer, size_t iouter, const field::BosOnv &inner, size_t iinner) {
        ASSERT_EQ(outer, chk_inds[iouter]);
        ASSERT_EQ(inner, chk_inds[iinner]);
    };
    mbf_foreach::bos::Pair<mbf_foreach::bos::GeneralOpen> foreach(hs);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(foreach.m_niter, chk_inds.size()*chk_inds.size());
}

TEST(MbfForeach, FrmBosGeneralOpen) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::general::chk_inds();
    const auto bos_chk_inds = bos::general_open::chk_inds();
    const sys::frm::Basis frm_hs(4, 3);
    const BosHilbertSpace bos_hs(0, 3, false, 2);

    size_t iiter = 0ul;
    auto fn = [&](const field::FrmBosOnv &field) {
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
        ++iiter;
    };
    mbf_foreach::frm::General outer(frm_hs);
    mbf_foreach::bos::GeneralOpen inner(bos_hs);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::General, mbf_foreach::bos::GeneralOpen> foreach(outer, inner);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, frm_chk_inds.size() * bos_chk_inds.size());
    ASSERT_EQ(foreach.m_niter, frm_chk_inds.size()*bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosSpins) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::spins::chk_inds();
    const auto bos_chk_inds = bos::general_open::chk_inds();
    const sys::frm::Basis frm_hs(4, 4, 0);
    const BosHilbertSpace bos_hs(0, 3, false, 2);

    size_t iiter = 0ul;
    auto fn = [&](const field::FrmBosOnv &field) {
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
        ++iiter;
    };
    mbf_foreach::frm::Spins outer(frm_hs);
    mbf_foreach::bos::GeneralOpen inner(bos_hs);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::Spins, mbf_foreach::bos::GeneralOpen> foreach(outer, inner);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, frm_chk_inds.size() * bos_chk_inds.size());
    ASSERT_EQ(foreach.m_niter, frm_chk_inds.size()*bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosSpinsEarlyExit) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::spins::chk_inds();
    const auto bos_chk_inds = bos::general_open::chk_inds();
    const sys::frm::Basis frm_hs(4, 4, 0);
    const BosHilbertSpace bos_hs(0, 3, false, 2);

    size_t iiter = 0ul;
    auto fn = [&](const field::FrmBosOnv &field) {
        if (iiter==32) throw ExitLoop();
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
        ++iiter;
    };
    mbf_foreach::frm::Spins outer(frm_hs);
    mbf_foreach::bos::GeneralOpen inner(bos_hs);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::Spins, mbf_foreach::bos::GeneralOpen> foreach(outer, inner);

    try { ASSERT_ANY_THROW(foreach.loop_fn(fn)); }
    catch (const ExitLoop &) {}
    ASSERT_EQ(iiter, 32);
    ASSERT_EQ(foreach.m_niter, frm_chk_inds.size()*bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosMs2Conserve) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::ms2_conserve::chk_inds();
    const auto bos_chk_inds = bos::general_open::chk_inds();
    const sys::frm::Basis frm_hs(5, 4, 1);
    const BosHilbertSpace bos_hs(0, 3, false, 2);

    size_t iiter = 0ul;
    auto fn = [&](const field::FrmBosOnv &field) {
        ASSERT_LT(iiter, frm_chk_inds.size() * bos_chk_inds.size());
        auto iiter_frm = iiter / bos_chk_inds.size();
        auto iiter_bos = iiter - iiter_frm * bos_chk_inds.size();
        ASSERT_TRUE(field.m_frm == frm_chk_inds[iiter_frm]);
        ASSERT_TRUE(field.m_bos == bos_chk_inds[iiter_bos]);
        ++iiter;
    };
    mbf_foreach::frm::Ms2Conserve outer(frm_hs);
    mbf_foreach::bos::GeneralOpen inner(bos_hs);
    mbf_foreach::frm_bos::Product<mbf_foreach::frm::Ms2Conserve, mbf_foreach::bos::GeneralOpen> foreach(outer, inner);
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    ASSERT_EQ(iiter, frm_chk_inds.size() * bos_chk_inds.size());
    ASSERT_EQ(foreach.m_niter, frm_chk_inds.size()*bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosMs2ConservePair) {
    using namespace mbf_foreach_test;
    const auto frm_chk_inds = frm::ms2_conserve::chk_inds();
    const auto bos_chk_inds = bos::general_open::chk_inds();
    const sys::frm::Basis frm_hs(5, 4, 1);
    const BosHilbertSpace bos_hs(0, 3, false, 2);

    auto fn = [&frm_chk_inds, &bos_chk_inds]
            (const field::FrmBosOnv &outer, size_t iouter, const field::FrmBosOnv &inner, size_t iinner) {
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
    typedef mbf_foreach::frm_bos::Product<mbf_foreach::frm::Ms2Conserve, mbf_foreach::bos::GeneralOpen> product_t;
    mbf_foreach::frm_bos::Pair<product_t> foreach({frm_hs, bos_hs});
    ASSERT_NO_THROW(foreach.loop_fn(fn));
    auto n = frm_chk_inds.size()*bos_chk_inds.size();
    ASSERT_EQ(foreach.m_niter, n*n);
}