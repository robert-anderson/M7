//
// Created by rja on 24/05/2020.
//

#include "InputOptions.h"

const std::string InputOptions::description =
        "\nM7: Many-body Stochastic Expectation Value Estimation Networks\n"
        "Command line interface\n";

InputOptions::InputOptions(CLI::App &app) : m_app(app) {

    add_option("-f,--fcidump_path", fcidump_path,
               "path to the FCIDUMP file, this is read-only");

    add_option("-O,--stats_path", fcidump_path,
               "path to the file to which MC cycle statistics will be output");

    add_flag("-E,--exact_propagation", exact_propagation,
             "perform fully deterministic projector FCI");

    add_option("-i,--nwalker_initial", nwalker_initial,
               "sum of walker magnitudes with which to initialize the populations.");

    add_option("-n,--nwalker_target", nwalker_target,
               "sum of walker magnitudes at which to begin varying the diagonal shift", true);

    add_option("-a,--nadd_initiator", nadd_initiator,
               "number of walkers defining the initiator threshold");

    add_option("-b,--max_bloom", max_bloom,
               "largest acceptable spawning bloom");

    add_option("-R,--prng_seed", prng_seed,
               "seed value for the mt19937 PRNG");

    add_option("--prng_ngen", prng_ngen,
               "number of PRNGs to batch generate with mt19937");

    add_option("-s,--ndet_semistoch", ndet_semistoch,
               "number of determinants selected to comprise the deterministic subspace");

    add_option("-z,--spin_restrict", spin_restrict,
               "difference in occupation of spin orbitals 0 and 1 in CI space for a spin-conserving hamiltonian");

    add_option("-W,--walker_factor_initial", walker_factor_initial,
               "number of rows initially allocated in the wavefunction store table as a multiple of the target walker number");

    add_option("-B,--buffer_factor_initial", buffer_factor_initial,
               "number of rows initially allocated in each segment of the wavefunction communicate buffer table as a multiple of the target walker number");

    add_option("-K,--nload_balance_block", nload_balance_block,
               "number of blocks per process to use for load balancing determinants among processes");

    add_option("-m,--min_spawn_mag", min_spawn_mag,
               "minimum spawn magnitude (stochastic threshold spawned weights about value)");

    add_option("-t,--tau_initial", tau_initial,
               "initial timestep");

    add_flag("--dynamic_tau", dynamic_tau,
             "update timestep dynamically");

    add_option("--nenough_spawns_for_dynamic_tau", nenough_spawns_for_dynamic_tau,
               "number of spawns logged for excitation type magnitudes to be used in tau update");

    add_option("-S,--shift_initial", shift_initial,
               "initial diagonal shift relative to the automatically assumed shift given by the reference energy.");

    add_option("-d,--shift_damp", shift_damp,
               "damping factor regulating the shift updates");

    add_option("-A,--shift_update_period", shift_update_period,
               "number of cycles between shift updates");

    add_option("-N,--ncycle", ncycle,
               "number of cycles to execute before exit.");
}