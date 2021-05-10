//
// Created by Robert John Anderson on 2020-02-08.
//

#include "Options.h"
#include "src/core/parallel/MPIAssert.h"

bool Options::init() {
    if (consts::float_is_zero(min_death_mag)) min_death_mag = min_spawn_mag;
    if (!defs::enable_bosons){
        MPI_REQUIRE(boson_coupling==0.0, "Boson coupling parameter is non-zero but bosons are compile time disabled");
        MPI_REQUIRE(boson_frequency==0.0, "Boson frequency parameter is non-zero but bosons are compile time disabled");
    }
    if (!defs::enable_mevs){
        MPI_REQUIRE(rdm_rank==0, "RDM rank is non-zero but MEVs are compile time disabled");
    }
    return true;
}
