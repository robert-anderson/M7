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
        auto exact_diag_energy = DenseHamiltonian(*fciqmc_calculation.m_ham).diagonalize().m_evals[0];
        ASSERT_FLOAT_EQ(consts::real(final_energy), exact_diag_energy);
    }
}

TEST(FciqmcCalculation, StochasticPropagation){
    Options options;
    options.fcidump_path = defs::assets_root+"/RHF_N2_6o6e/FCIDUMP";
    options.tau_initial = 0.05;
    options.nwalker_target = 150000;
    options.nload_balance_block = 5;
    options.ncycle = 4000;
    options.prng_seed = 12;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
        auto energy_mean_std = stat_utils::quotient(num, den);
        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
        ASSERT_TRUE(consts::floats_nearly_equal(
            consts::real(energy_mean_std.first), -108.8113865756313, 1e-3));
    }
}

TEST(FciqmcCalculation, StochasticPropagation4c){
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();
    Options options;
    options.fcidump_path = defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP";
    options.tau_initial = 0.05;
    options.nwalker_target = 150000;
    options.nload_balance_block = 5;
    options.ncycle = 4000;
    options.prng_seed = 12;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
        auto energy_mean_std = stat_utils::quotient(num, den);
        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
        ASSERT_TRUE(consts::floats_nearly_equal(
                consts::real(energy_mean_std.first), -14.28882489, 1e-3));
    }
}


TEST(FciqmcCalculation, StochasticPropagation60orbs){
    Options options;
    options.fcidump_path = defs::assets_root+"/RHF_N2_CCPVTZ/FCIDUMP";
    options.tau_initial = 0.002;
    options.nwalker_target = 150000;
    options.ncycle = 4000;
    options.prng_seed = 12;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
        auto energy_mean_std = stat_utils::quotient(num, den);
        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
        //ASSERT_TRUE(consts::floats_nearly_equal(energy_mean_std.first, -108.8113865756313, 1e-3));
    }
}

TEST(FciqmcCalculation, SemiStochasticPropagation){
    Options options;
    options.fcidump_path = defs::assets_root+"/RHF_N2_6o6e/FCIDUMP";
    options.tau_initial = 0.05;
    options.prng_seed = 13;
    options.nwalker_target = 100000;
    options.ncycle = 100;
    options.do_semistochastic = true;
    options.niter_init_detsub = 1000;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
        auto energy_mean_std = stat_utils::quotient(num, den);
        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
        ASSERT_TRUE(consts::floats_nearly_equal(
            consts::real(energy_mean_std.first), -108.8113865756313, 1e-3));
    }
}