//
// Created by Robert J. Anderson on 24/05/2020.
//

#include <M7_lib/dynamics/StochasticPropagator.h>
#include <M7_lib/linalg/DenseHamiltonian.h>
#include "gtest/gtest.h"
#include "M7_lib/dynamics/FciqmcCalculation.h"

#if 0
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
    options.nload_balance_block_per_rank = 5;
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
    options.nload_balance_block_per_rank = 5;
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
    options.tau_initial = 0.24923054538193248;
    options.prng_seed = 13;
    options.nwalker_target = 100000;
    options.nwalker_initial = 100;
    options.ncycle = 10000;
    //options.exact_propagation = 1;
    options.do_semistochastic = true;
    options.ncycle_wait_detsub = 1000;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
//    if (mpi::i_am_root()) {
//        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
//        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
//        auto energy_mean_std = stat_utils::quotient(num, den);
//        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
//        ASSERT_TRUE(consts::floats_nearly_equal(
//            consts::real(energy_mean_std.first), -108.8113865756313, 1e-3));
//    }
}


TEST(FciqmcCalculation, SemiStochasticPropagation4Fold){
    Options options;
    options.fcidump_path = defs::assets_root+"/HF_DIRAC_4fold/FCIDUMP";
    options.shift_initial = 0.2;
    options.fcidump_spin_major = 1;
    options.tau_initial = 0.1;
    options.static_tau = true;
    options.prng_seed = 13;
    options.nwalker_target = 100000;
    options.nwalker_initial = 100;
    options.ncycle = 50000;
    options.do_semistochastic = true;
    options.nadd_thresh_semistoch = 1;
    options.ncycle_wait_detsub = 1000;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();
    if (mpi::i_am_root()) {
        auto num = fciqmc_calculation.m_stats_file->m_ref_proj_energy_num.mean_std(options.ncycle / 2);
        auto den = fciqmc_calculation.m_stats_file->m_ref_weight.mean_std(options.ncycle / 2);
        auto energy_mean_std = stat_utils::quotient(num, den);
        std::cout << std::setprecision(10) << energy_mean_std.first << std::endl;
    }
}



TEST(FciqmcCalculation, StochasticPropagation4cLarge){
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();
    Options options;
    options.fcidump_path = defs::assets_root+"/DHF_Kr_CCPVDZ/FCIDUMP";
    options.tau_initial = 0.01;
    options.nwalker_target = 1000000;
    options.nload_balance_block_per_rank = 1000;
    options.ncycle = 4000;
    options.walker_buffer_size_factor_initial = 2.0;
    options.spawn_buffer_size_factor_initial = 3.0;
    options.do_semistochastic = true;
    options.ncycle_wait_detsub = 3000;
    FciqmcCalculation fciqmc_calculation(options);
    fciqmc_calculation.execute();

}

#endif