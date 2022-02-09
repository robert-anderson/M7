//
// Created by anderson on 09/02/2022.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/basis/MbfForeach.h>
#include "gtest/gtest.h"

TEST(MbfForeach, FrmGeneral) {

#if 0
    const size_t nsite=4, nelec=4;
    const size_t nmode=3, nboson_max=2;
    buffered::FrmBosOnv mbf({nsite, nmode});

//    struct Foreach : mbf_foreach::FrmOnvGeneral {
//        Foreach(field::FrmOnv& mbf, size_t nelec): mbf_foreach::FrmOnvGeneral(mbf, nelec){}
//
//        void mbf_body(const field::FrmOnv &mbf) override {
//            std::cout << mbf << std::endl;
//        }
//    };

    struct Foreach : mbf_foreach::FrmBosOnvGeneral {
        Foreach(field::FrmBosOnv& mbf, size_t nelec, size_t nboson_max):
            mbf_foreach::FrmBosOnvGeneral(mbf, nelec, nboson_max){}

        void mbf_body(const field::FrmBosOnv &mbf) override {
            std::cout << m_mbf << std::endl;
        }
    };

    Foreach foreach(mbf, nelec, nboson_max);
    foreach.loop();
#endif
}