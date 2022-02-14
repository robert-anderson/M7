//
// Created by anderson on 09/02/2022.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/basis/MbfForeach.h>
#include "gtest/gtest.h"

namespace mbf_foreach_test {

    namespace frm {

        struct General : mbf_foreach::frm::General {
            const std::vector<defs::inds> m_chk_inds;

            General(field::FrmOnv *mbf = nullptr) : mbf_foreach::frm::General(3, 4, mbf), m_chk_inds(chk_inds()) {}
            General(const General &other, field::FrmOnv *mbf = nullptr) : General(mbf) {}

            void body() override {
                ASSERT_LT(iiter(), m_chk_inds.size());
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }

            static std::vector<defs::inds> chk_inds() {
                return {
                    {0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}, {1, 2, 3, 4}, {0, 1, 2, 5},
                    {0, 1, 3, 5}, {0, 2, 3, 5}, {1, 2, 3, 5}, {0, 1, 4, 5}, {0, 2, 4, 5}, {1, 2, 4, 5},
                    {0, 3, 4, 5}, {1, 3, 4, 5}, {2, 3, 4, 5}
                };
            };
        };

        struct Spins : mbf_foreach::frm::Spins {
            const std::vector<defs::inds> m_chk_inds;

            Spins(field::FrmOnv *mbf = nullptr) : mbf_foreach::frm::Spins(4, 0, mbf), m_chk_inds(chk_inds()) {}
            Spins(const Spins &other, field::FrmOnv *mbf = nullptr) : Spins(mbf) {}

            void body() override {
                ASSERT_LT(iiter(), m_chk_inds.size());
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }
            
            static std::vector<defs::inds> chk_inds() {
                return {
                    {0, 1, 6, 7}, {0, 2, 5, 7}, {1, 2, 4, 7}, {0, 3, 5, 6}, {1, 3, 4, 6}, {2, 3, 4, 5}
                };
            }
        };

        struct Ms2Conserve : mbf_foreach::frm::Ms2Conserve {
            const std::vector<defs::inds> m_chk_inds;

            Ms2Conserve(field::FrmOnv *mbf = nullptr) :
                    mbf_foreach::frm::Ms2Conserve(4, 5, 1, mbf), m_chk_inds(chk_inds()) {}

            Ms2Conserve(const Ms2Conserve &other, field::FrmOnv *mbf = nullptr) : Ms2Conserve(mbf) {}

            void body() override {
                ASSERT_LT(iiter(), m_chk_inds.size());
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }
            
            static std::vector<defs::inds> chk_inds() {
                return {
                    {0, 1, 2, 4, 5}, {0, 1, 2, 4, 6}, {0, 1, 2, 5, 6}, {0, 1, 2, 4, 7}, {0, 1, 2, 5, 7}, {0, 1, 2, 6, 7},
                    {0, 1, 3, 4, 5}, {0, 1, 3, 4, 6}, {0, 1, 3, 5, 6}, {0, 1, 3, 4, 7}, {0, 1, 3, 5, 7}, {0, 1, 3, 6, 7},
                    {0, 2, 3, 4, 5}, {0, 2, 3, 4, 6}, {0, 2, 3, 5, 6}, {0, 2, 3, 4, 7}, {0, 2, 3, 5, 7}, {0, 2, 3, 6, 7},
                    {1, 2, 3, 4, 5}, {1, 2, 3, 4, 6}, {1, 2, 3, 5, 6}, {1, 2, 3, 4, 7}, {1, 2, 3, 5, 7}, {1, 2, 3, 6, 7}
                };
            }
        };
    }

    namespace bos {

        struct General : mbf_foreach::bos::General {
            const std::vector<defs::inds> m_chk_inds;

            General(field::BosOnv *mbf = nullptr) : mbf_foreach::bos::General(3, 2, mbf),
                m_chk_inds(chk_inds()) {}

            General(const General &other, field::BosOnv *mbf = nullptr) : General(mbf) {}

            void body() override {
                ASSERT_LT(iiter(), m_chk_inds.size());
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }            
            
            static std::vector<defs::inds> chk_inds() {
                return {
                    {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 0}, {0, 1, 1}, {0, 1, 2}, {0, 2, 0}, {0, 2, 1}, {0, 2, 2},
                    {1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {1, 2, 0}, {1, 2, 1}, {1, 2, 2},
                    {2, 0, 0}, {2, 0, 1}, {2, 0, 2}, {2, 1, 0}, {2, 1, 1}, {2, 1, 2}, {2, 2, 0}, {2, 2, 1}, {2, 2, 2}
                };
            }
        };
    }

    namespace frm_bos {

        template<typename frm_t, typename bos_t>
        struct Product : mbf_foreach::frm_bos::Product<frm_t, bos_t> {
            typedef mbf_foreach::frm_bos::Product<frm_t, bos_t> base_t;
            const std::vector<defs::inds> m_frm_chk_inds;
            const std::vector<defs::inds> m_bos_chk_inds;

            Product() : base_t({}, {}),
                m_frm_chk_inds(frm_t::chk_inds()), m_bos_chk_inds(bos_t::chk_inds()){}

            using base_t::iiter;
            using base_t::value;
            void body() override {
                ASSERT_LT(iiter(), m_frm_chk_inds.size()*m_bos_chk_inds.size());
                auto iiter_frm = iiter()/m_bos_chk_inds.size();
                auto iiter_bos = iiter()-iiter_frm*m_bos_chk_inds.size();
                ASSERT_TRUE(value().m_frm == m_frm_chk_inds[iiter_frm]);
                ASSERT_TRUE(value().m_bos == m_bos_chk_inds[iiter_bos]);
            }
        };
    }
}

TEST(MbfForeach, FrmGeneral) {
    using namespace mbf_foreach_test::frm;
    General foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_chk_inds.size());
}

TEST(MbfForeach, FrmSpins) {
    using namespace mbf_foreach_test::frm;
    Spins foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_chk_inds.size());
}

TEST(MbfForeach, FrmMs2Conserve) {
    using namespace mbf_foreach_test::frm;
    Ms2Conserve foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_chk_inds.size());
}

TEST(MbfForeach, BosGeneral) {
    using namespace mbf_foreach_test::bos;
    General foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_chk_inds.size());
}

TEST(MbfForeach, FrmBosGeneral) {
    using namespace mbf_foreach_test;
    frm_bos::Product<frm::General, bos::General> foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_frm_chk_inds.size()*foreach.m_bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosSpins) {
    using namespace mbf_foreach_test;
    frm_bos::Product<frm::Spins, bos::General> foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_frm_chk_inds.size()*foreach.m_bos_chk_inds.size());
}

TEST(MbfForeach, FrmBosMs2Conserve) {
    using namespace mbf_foreach_test;
    frm_bos::Product<frm::Ms2Conserve, bos::General> foreach;
    foreach.loop();
    ASSERT_EQ(foreach.iiter()+1, foreach.m_frm_chk_inds.size()*foreach.m_bos_chk_inds.size());
}