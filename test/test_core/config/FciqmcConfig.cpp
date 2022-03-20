//
// Created by rja on 25/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/config/FciqmcConfig.h"

TEST(FciqmcConfig, Test){
    fciqmc_config::Document doc(nullptr);
    std::cout <<
              doc.help_string()
            << std::endl;
}