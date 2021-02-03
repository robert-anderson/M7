//
// Created by Robert John Anderson on 2020-02-08.
//

#include "Options.h"

bool Options::init() {
    if (consts::float_is_zero(min_death_mag)) min_death_mag = min_spawn_mag;
    return true;
}
