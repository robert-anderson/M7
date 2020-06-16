//
// Created by rja on 24/05/2020.
//

#include <src/core/dynamics/StochasticPropagator.h>
#include <src/core/linalg/DenseHamiltonian.h>
#include <src/core/linalg/EigenSolver.h>
#include "gtest/gtest.h"
#include "src/core/dynamics/FciqmcCalculation.h"

TEST(FciqmcCalculation, ExactPropagation){
    Options options;
    options.fcidump_path = defs::assets_root+"/RHF_N2_6o6e/FCIDUMP";
    options.tau_initial = 0.05;
    options.nwalker_target = 1000;
    options.ncycle = 1000;
    options.exact_propagation = 1;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto final_energy = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.m_series.back() /
                            fciqmc_calculation.m_stats_file->m_ref_weight.m_series.back();
        // check exact propagation energy against exact eigensolver
        auto exact_diag_energy = DenseHamiltonian(*fciqmc_calculation.m_ham).diagonalize().m_evals(0);
        ASSERT_FLOAT_EQ(consts::real(final_energy), exact_diag_energy);
    }
}

TEST(FciqmcCalculation, StochasticPropagation){
    Options options;
    options.fcidump_path = defs::assets_root+"/RHF_N2_6o6e/FCIDUMP";
    options.tau_initial = 0.05;
    options.nwalker_target = 10000;
    options.ncycle = 20000;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
        auto energy_mean_std = stat_utils::quotient(num, den);
        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
        ASSERT_TRUE(consts::floats_nearly_equal(energy_mean_std.first, -108.8113865756313, 1e-3));
    }
}