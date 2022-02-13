//
// Created by anderson on 09/02/2022.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/basis/MbfForeach.h>
#include <src/core/table/BufferedTable.h>
#include "gtest/gtest.h"

namespace mbf_foreach_test {
    namespace spins {
        struct Foreach : mbf_foreach::FrmOnvSpins {
            std::vector<defs::inds> m_chk_inds;

            Foreach() : mbf_foreach::FrmOnvSpins(4, 0) {
                m_chk_inds.reserve(6);
                m_chk_inds.push_back({0, 1, 6, 7});
                m_chk_inds.push_back({0, 2, 5, 7});
                m_chk_inds.push_back({1, 2, 4, 7});
                m_chk_inds.push_back({0, 3, 5, 6});
                m_chk_inds.push_back({1, 3, 4, 6});
                m_chk_inds.push_back({2, 3, 4, 5});
            }

            void body() override {
                std::cout << m_mbf << std::endl;
                std::cout << iiter() << std::endl;
                std::cout << m_chk_inds[iiter()] << std::endl;
                ASSERT_TRUE(m_mbf == m_chk_inds[iiter()]);
            }
        };
    }
}

TEST(MbfForeach, FrmOnvGeneral) {

    struct Foreach : mbf_foreach::FrmOnvGeneral {
        Foreach(size_t nsite, size_t nelec):
                mbf_foreach::FrmOnvGeneral(nsite, nelec){}

        void body() override {
            std::cout << m_mbf << std::endl;
        }
    };

    Foreach foreach(6, 6);
    foreach.loop();
}

TEST(MbfForeach, FrmOnvSpins) {
    using namespace mbf_foreach_test::spins;
    Foreach foreach;
    foreach.loop();
}


TEST(MbfForeach, FrmOnvSzConserve) {

    struct Foreach : mbf_foreach::FrmOnvSzConserve {
        Foreach(size_t nsite, size_t nelec):
                mbf_foreach::FrmOnvSzConserve(nsite, nelec / 2){}

        void body() override {
            std::cout << m_mbf << std::endl;
        }
    };

    Foreach foreach(6, 6);
    foreach.loop();
}






