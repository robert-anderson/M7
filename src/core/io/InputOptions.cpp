//
// Created by rja on 24/05/2020.
//

#include "InputOptions.h"

const std::string InputOptions::program_description =
        "\nM7: Many-body Stochastic Expectation Value Estimation Networks\n"
        "Command line interface\n";

InputOptions::InputOptions(CLI::App &app) : m_app(app) {

    add_option("-f,--fcidump_path", fcidump_path,
               "path to the FCIDUMP file, this is read-only");

    add_option("-D,--initial-reference-det", initial_reference_det,
               "initial reference determinant");

    add_option("-y, --reference_redefinition_thresh", reference_redefinition_thresh,
               "the reference will change to any determinant whose weight exceeds the reference weight by this factor");

    add_flag("--spin_major", fcidump_spin_major,
               "if true, spin-resolved FCIDUMP orders the spin orbitals aaa...bbb..., and ababab... if false.");

    add_option("-O,--stats_path", stats_path,
               "path to the file to which MC cycle statistics will be output");

    add_flag("-E,--exact_propagation", exact_propagation,
             "perform fully deterministic projector FCI");

    add_set("-e,--excit_gen", excit_gen, {"pchb", "uniform"},
             "excitation generator to use for stochastic connections");

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

    add_option("--walker_fraction_semistoch", walker_fraction_semistoch,
               "determinant walker occupation as a fraction of total N_W for inclusion in deterministic subspace");

    add_option("--nadd_thresh_semistoch", nadd_thresh_semistoch,
               "determinant walker occupation in units of nadd_initiator for inclusion in deterministic subspace");

    add_option("-z,--spin_restrict", spin_restrict,
               "difference in occupation of spin orbitals 0 and 1 in CI space for a spin-conserving hamiltonian");

    add_option("-W,--walker_factor_initial", walker_factor_initial,
               "number of rows initially allocated in the wavefunction store table as a multiple of the target walker number");

    add_option("-B,--buffer_factor_initial", buffer_factor_initial,
               "number of rows initially allocated in each segment of the wavefunction communicate buffer table as a multiple of the target walker number");

    add_option("-K,--nload_balance_block_per_rank", nload_balance_block_per_rank,
               "number of blocks per process to use for load balancing determinants among processes");

    add_option("-m,--min_spawn_mag", min_spawn_mag,
               "minimum spawn magnitude (stochastic threshold spawned weights about value)");

    add_option("-t,--tau_initial", tau_initial,
               "initial timestep");

    add_flag("--static_tau", static_tau,
             "disable dynamic timestep update");

    add_option("--min_excit_class_prob", min_excit_class_prob,
             "prevent the probability of drawing an excitatation class falling below this threshold");

    add_option("--nenough_spawns_for_dynamic_tau", nenough_spawns_for_dynamic_tau,
               "number of spawns logged for excitation type magnitudes to be used in tau update");

    add_option("--shift_initial", shift_initial,
               "initial diagonal shift relative to the automatically assumed shift given by the reference energy.");

    add_option("-d,--shift_damp", shift_damp,
               "damping factor regulating the shift updates");

    add_option("-A,--shift_update_period", shift_update_period,
               "number of cycles between shift updates");

    add_option("-N,--ncycle", ncycle,
               "number of cycles to execute before exit");

    add_flag("-S,--semistochastic", do_semistochastic,
             "enable semistochastic adaptation");

    add_option("--ncycle_init_detsub", ncycle_init_detsub,
               "number of cycles after start of variable shift epoch begin semistochastic epoch");

    add_flag("--calc_mk_walker_sums", calc_mk_walker_sums,
             "accumulate and output average walker occupations of Kramers sectors");

}