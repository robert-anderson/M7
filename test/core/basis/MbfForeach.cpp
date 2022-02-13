//
// Created by anderson on 09/02/2022.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/basis/MbfForeach.h>
#include "gtest/gtest.h"

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

    struct Foreach : mbf_foreach::FrmOnvSpins {
        Foreach(size_t nsite):
                mbf_foreach::FrmOnvSpins(nsite, 0){}

        void body() override {
            std::cout << m_mbf << std::endl;
        }
    };

    Foreach foreach(6);
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






