//
// Created by Robert John Anderson on 2020-02-08.
//

#ifndef M7_OPTIONS_H
#define M7_OPTIONS_H

#include <string>

struct Options {
    std::string fcidump_path = "FCIDUMP";
    std::string initial_reference_det;
    double reference_redefinition_thresh = 2.0;
    bool fcidump_spin_major = false;
    std::string stats_path = "M7.stats";
    std::string parallel_stats_path = "M7.rank";
    bool exact_propagation = false;
    std::string excit_gen = "pchb";
    double nwalker_initial = 1.0;
    double nwalker_target = 0.0;
    double nadd_initiator = 3.0;
    double max_bloom = 3.0;
    size_t prng_seed = 12;
    size_t prng_ngen = 10000;
    size_t ndet_semistoch = 0;
    double walker_fraction_semistoch = 1.0;
    double nadd_thresh_semistoch = 0.0;
    size_t spin_restrict = 0;
    double walker_factor_initial = 1.0;
    double buffer_factor_initial = 10.0;
    double buffer_expansion_factor=0.5;
    double min_spawn_mag = 0.0;
    size_t nload_balance_block_per_rank = 10;
    size_t load_balance_period = 10;
    double tau_initial = 0.05;
    bool static_tau = false;
    double min_excit_class_prob = 1e-3;
    size_t nenough_spawns_for_dynamic_tau = 1000;
    double shift_initial = 0.0;
    double shift_damp = 0.05;
    size_t shift_update_period = 1;
    size_t ncycle = ~0ul;
    bool do_semistochastic = false;
    size_t ncycle_init_detsub = 1000;
    bool calc_mk_walker_sums = false;
    size_t nboson_max = 0ul;
    double boson_coupling = 0.0;
    double boson_frequency = 0.0;

    bool validate() const;

};


#endif //M7_OPTIONS_H
