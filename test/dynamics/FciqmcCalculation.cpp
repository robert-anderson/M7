//
// Created by rja on 24/05/2020.
//

#include "gtest/gtest.h"
#include "src/core/dynamics/FciqmcCalculation.h"

TEST(FciqmcCalculation, DHF_Be_STO3G){
    Options options;
    options.fcidump_path = defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP";
    options.nwalker_target = 1000;
    //options.ncycle = 1;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
}