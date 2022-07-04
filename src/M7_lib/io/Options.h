//
// Created by Robert John Anderson on 2020-02-08.
//

#ifndef M7_OPTIONS_H
#define M7_OPTIONS_H

#include <string>
#include <M7_lib/defs.h>

#if 0
struct Options {
    str_t fcidump_path = "FCIDUMP";
    str_t initial_reference_det;
    double reference_redefinition_thresh = 10.0;
    bool fcidump_spin_major = false;
    str_t stats_path = "M7.stats";
    bool parallel_stats = false;
    bool exact_propagation = false;
    str_t excit_gen = "pchb";
    bool spf_uniform_twf = false;
    bool spf_weighted_twf = false;
    double nwalker_initial = 1.0;
    double nwalker_target = 0.0;
    double nadd_initiator = 3.0;
    double max_bloom = 3.0;
    uint_t nroot = 1;
    uint_t prng_seed = 12;
    uint_t prng_ngen = 10000;
    uint_t ndet_semistoch = 0;
    double walker_fraction_semistoch = 1.0;
    double nadd_thresh_semistoch = 0.0;
    uint_t spin_restrict = 0;
    double walker_buffer_size_factor_initial = 1.0;
    double walker_buffer_expansion_factor = 0.5;
    double mev_buffer_expansion_factor = 0.5;
    double spawn_buffer_size_factor_initial = 10.0;
    double min_spawn_mag = 0.4;
    double min_death_mag = 0.0;
    uint_t nload_balance_block_per_rank = 20;
    uint_t load_balance_period = 10;
    double acceptable_load_imbalance = 0.05;
    double tau_initial = 0.05;
    bool static_tau = false;
    double min_excit_class_prob = 1e-3;
    uint_t nenough_spawns_for_dynamic_tau = 1000;
    double shift_initial = 0.0;
    double shift_damp = 0.05;
    uint_t shift_update_period = 1;
    uint_t ncycle_reweight_lookback = 0;
    uint_t ncycle_wait_reweight = 2000ul;
    uint_t ncycle = ~0ul;
    uint_t ncycle_wait_mevs = 1000ul;
    uint_t ncycle_accumulate_mevs = ~0ul;
    uint_t ncycle_mev_period = 50ul;
    uint_t ncycle_shift_average_period = 1000ul;
    uint_t max_rank_average_coeff = 0;
    bool output_mevs_periodically = false;
    bool mev_mixed_estimator = false;
    bool consolidate_spawns = false;
    bool explicit_hf_conn_mevs = true;
    bool do_semistochastic = false;
    uint_t ncycle_wait_detsub = 100;
    bool calc_mk_walker_sums = false;
    uint_t nboson_max = 0ul;
    double boson_coupling = 0.0;
    double boson_frequency = 0.0;
    double spf_twf_fermion_factor = 1.0;
    double spf_twf_boson_factor = 1.0;
    prob_t psingle_initial = 0.0;
    uint_t rdm_rank = 0;
    str_t write_hdf5_fname = "";
    str_t read_hdf5_fname = "";
    bool replicate = false;

    bool init();

};
#endif


#endif //M7_OPTIONS_H
