//
// Created by Robert John Anderson on 2020-02-08.
//

#include "InputOptions.h"

const std::string InputOptions::description =
    "\nM7: Many-body Stochastic Expectation Value Estimation Networks\n"
    "Command line interface\n";

InputOptions::InputOptions(CLI::App &app) : m_app(app) {

    add_option("-f,--fcidump_path", fcidump_path,
               "Path to the FCIDUMP file, this is read-only");

    add_option("-O,--stats_path", fcidump_path,
               "Path to the file to which MC cycle statistics will be output");

    add_option("-i,--nwalker_initial", nwalker_initial,
               "sum of walker magnitudes with which to initialize the populations.");

    add_option("-n,--nwalker_target", nwalker_target,
               "sum of walker magnitudes at which to begin varying the diagonal shift");

    add_option("-a,--nadd_initiator", nadd_initiator,
               "Number of walkers defining the initiator threshold");

    add_option("-R,--prng_seed", prng_seed,
               "Seed value for the mt19937 PRNG");

    add_option("-s,--ndet_semistoch", ndet_semistoch,
               "Number of determinants selected to comprise the deterministic subspace");

    add_option("-z,--spin_restrict", spin_restrict,
               "Difference in occupation of spin orbitals 0 and 1 in CI space for a spin-conserving hamiltonian");

    add_option("-W,--walker_factor_initial", walker_factor_initial,
               "number of rows initially allocated in the wavefunction store table as a multiple of the target walker number");

    add_option("-B,--buffer_factor_initial", buffer_factor_initial,
               "number of rows initially allocated in each segment of the wavefunction send buffer table as a multiple of the target walker number");

    add_option("-K,--nload_balance_block", nload_balance_block,
               "number of blocks per process to use for load balancing determinants among processes");

    add_option("-t,--tau_initial", tau_initial,
            "initial timestep");

    add_option("-S,--shift_initial", shift_initial,
            "initial diagonal shift relative to the automatically assumed shift given by the reference energy.");

    add_option("-d,--shift_damp", shift_damp,
               "damping factor regulating the shift updates");

    add_option("-A,--shift_update_period", shift_update_period,
               "number of cycles between shift updates");

    add_option("-N,--ncycle", ncycle,
               "number of cycles to execute before exit.");



}

bool InputOptions::validate() const {
    return true;
}
