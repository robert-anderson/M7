//
// Created by rja on 25/06/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/config/FciqmcConfig.h"

TEST(FciqmcConfig, Test){
    conf::Document doc(nullptr);
    std::cout <<
              doc.help_string()
            << std::endl;
}