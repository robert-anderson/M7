//
// Created by Robert John Anderson on 2020-02-08.
//

#include "InputOptions.h"

const std::string InputOptions::description =
    "\nM7: Many-body Stochastic Expectation Value Estimation Networks\n"
    "Command line interface\n";

InputOptions::InputOptions(CLI::App &app) : m_app(app) {

    add_option("-f,--fcidump-path", fcidump_path,
               "Path to the FCIDUMP file, this is read-only");

    add_option("-i,--nwalker-initial", nwalker_initial,
               "Number of walkers with which the populations will be initialized");

    add_option("-n,--nwalker-target", nwalker_target,
               "Number of walkers at which the population will be stabilised after the growth phase", true);

    add_option("-a,--nadd-initiator", nadd_initiator,
               "Number of walkers defining the initiator threshold");

    add_option("-s,--ndet-semistoch", ndet_semistoch,
               "Number of determinants selected to comprise the deterministic subspace");

    add_option("-o,--spin-level", spin_level,
               "n/2 (even nelec), or (n-1)/2 (odd nelec) where n is number of open shells in the reference determinant");

    add_option("-S,--store_factor_initial", store_factor_initial,
               "number of rows initially allocated in the wavefunction store table as a multiple of the target walker number");

    add_option("-B,--buffer_factor_initial", buffer_factor_initial,
               "number of rows initially allocated in each segment of the wavefunction send buffer table as a multiple of the target walker number");

    add_option("-T,--buffer_temp_factor_initial", buffer_temp_factor_initial,
               "number of rows allocated in each segment of the temporary wavefunction send buffer table as a multiple of the target walker number");

    add_option("-K,--nload_balance_block", nload_balance_block,
               "number of blocks per process to use for load balancing determinants among processes");


}

bool InputOptions::validate() const {
    return true;
}
