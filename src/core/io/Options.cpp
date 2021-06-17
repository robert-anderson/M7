//
// Created by Robert John Anderson on 2020-02-08.
//

#include "Options.h"
#include "src/core/parallel/MPIAssert.h"

bool Options::init() {
    /*
     * Sanity checks on input options, and setting defaults of options whose defaults depend on the value of other
     * options.
     */
    if (consts::float_is_zero(min_death_mag)) min_death_mag = min_spawn_mag;

    if (ncycle_wait_reweight < ncycle_shift_average_period) {
        log::warn("ncycle_wait_reweight cannot be less than ncycle_shift_average_period. Setting them equal.");
        ncycle_wait_reweight = ncycle_shift_average_period;
    }

    if (!defs::enable_bosons){
        REQUIRE_EQ_ALL(boson_coupling, 0.0, "Boson coupling parameter is non-zero but bosons are compile time disabled");
        REQUIRE_EQ_ALL(boson_frequency, 0.0, "Boson frequency parameter is non-zero but bosons are compile time disabled");
    }
    if (!defs::enable_mevs){
        REQUIRE_EQ_ALL(rdm_rank, 0, "RDM rank is non-zero but MEVs are compile time disabled");
    }
    return true;
}
