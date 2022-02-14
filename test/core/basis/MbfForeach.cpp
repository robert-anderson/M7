//
// Created by anderson on 09/02/2022.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/basis/MbfForeach.h>
#include "gtest/gtest.h"

namespace mbf_foreach_test {

    namespace frm {

        struct General : mbf_foreach::frm::General {
            std::vector<defs::inds> m_chk_inds;

            General() : mbf_foreach::frm::General(3, 4) {
                m_chk_inds = {
                        {0, 1, 2, 3},
                        {0, 1, 2, 4},
                        {0, 1, 3, 4},
                        {0, 2, 3, 4},
                        {1, 2, 3, 4},
                        {0, 1, 2, 5},
                        {0, 1, 3, 5},
                        {0, 2, 3, 5},
                        {1, 2, 3, 5},
                        {0, 1, 4, 5},
                        {0, 2, 4, 5},
                        {1, 2, 4, 5},
                        {0, 3, 4, 5},
                        {1, 3, 4, 5},
                        {2, 3, 4, 5}
                };
            }

            void body() override {
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }
        };

        struct Spins : mbf_foreach::frm::Spins {
            std::vector<defs::inds> m_chk_inds;

            Spins() : mbf_foreach::frm::Spins(4, 0) {
                m_chk_inds = {
                        {0, 1, 6, 7},
                        {0, 2, 5, 7},
                        {1, 2, 4, 7},
                        {0, 3, 5, 6},
                        {1, 3, 4, 6},
                        {2, 3, 4, 5}
                };
            }

            void body() override {
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }
        };

        struct Ms2Conserve : mbf_foreach::frm::Ms2Conserve {
            std::vector<defs::inds> m_chk_inds;

            Ms2Conserve() : mbf_foreach::frm::Ms2Conserve(4, 5, 1) {
                m_chk_inds = {
                        {0, 1, 2, 4, 5},
                        {0, 1, 2, 4, 6},
                        {0, 1, 2, 5, 6},
                        {0, 1, 2, 4, 7},
                        {0, 1, 2, 5, 7},
                        {0, 1, 2, 6, 7},

                        {0, 1, 3, 4, 5},
                        {0, 1, 3, 4, 6},
                        {0, 1, 3, 5, 6},
                        {0, 1, 3, 4, 7},
                        {0, 1, 3, 5, 7},
                        {0, 1, 3, 6, 7},

                        {0, 2, 3, 4, 5},
                        {0, 2, 3, 4, 6},
                        {0, 2, 3, 5, 6},
                        {0, 2, 3, 4, 7},
                        {0, 2, 3, 5, 7},
                        {0, 2, 3, 6, 7},

                        {1, 2, 3, 4, 5},
                        {1, 2, 3, 4, 6},
                        {1, 2, 3, 5, 6},
                        {1, 2, 3, 4, 7},
                        {1, 2, 3, 5, 7},
                        {1, 2, 3, 6, 7}
                };
            }

            void body() override {
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }
        };
    }

    namespace bos {

        struct General : mbf_foreach::bos::General {
            std::vector<defs::inds> m_chk_inds;

            General() : mbf_foreach::bos::General(3, 2) {
                m_chk_inds = {
                        {0, 0, 0},
                        {0, 0, 1},
                        {0, 0, 2},
                        {0, 1, 0},
                        {0, 1, 1},
                        {0, 1, 2},
                        {0, 2, 0},
                        {0, 2, 1},
                        {0, 2, 2},

                        {1, 0, 0},
                        {1, 0, 1},
                        {1, 0, 2},
                        {1, 1, 0},
                        {1, 1, 1},
                        {1, 1, 2},
                        {1, 2, 0},
                        {1, 2, 1},
                        {1, 2, 2},

                        {2, 0, 0},
                        {2, 0, 1},
                        {2, 0, 2},
                        {2, 1, 0},
                        {2, 1, 1},
                        {2, 1, 2},
                        {2, 2, 0},
                        {2, 2, 1},
                        {2, 2, 2}
                };
            }

            void body() override {
                ASSERT_TRUE(value() == m_chk_inds[iiter()]);
            }
        };
    }

    namespace frm_bos {
        typedef mbf_foreach::frm_bos::Product<frm::Ms2Conserve, bos::General> product_base_t;
        struct ProductMs2Conserve : product_base_t {
            std::vector<defs::inds> m_chk_inds;
            ProductMs2Conserve() : product_base_t(frm::Ms2Conserve(), bos::General()) {}
        };
    }
}

TEST(MbfForeach, FrmGeneral) {
    using namespace mbf_foreach_test::frm;
    General().loop();
}

TEST(MbfForeach, FrmSpins) {
    using namespace mbf_foreach_test::frm;
    Spins().loop();
}

TEST(MbfForeach, FrmMs2Conserve) {
    using namespace mbf_foreach_test::frm;
    Ms2Conserve().loop();
}

TEST(MbfForeach, BosGeneral) {
    using namespace mbf_foreach_test::bos;
    General().loop();
}

TEST(MbfForeach, FrmBosProduct) {
    using namespace mbf_foreach_test::frm_bos;
    ProductMs2Conserve().loop();
}