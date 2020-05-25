//
// Created by rja on 24/05/2020.
//

#include <src/core/dynamics/StochasticPropagator.h>
#include "gtest/gtest.h"
#include "src/core/dynamics/FciqmcCalculation.h"

TEST(FciqmcCalculation, DHF_Be_STO3G){
    Options options;
    options.fcidump_path = defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP";
    options.tau_initial = 0.05;
    options.nwalker_target = 1000;
    options.ncycle = 10;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    utils::print(fciqmc_calculation.m_stats_file.m_ref_proj_energy_num.m_series);
}