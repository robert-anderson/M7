//
// Created by rja on 27/02/2020.
//

#include <src/defs.h>
#include "SemistochasticPropagator.h"

SemistochasticPropagator::SemistochasticPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                                                   const RankAllocator<Determinant> &rank_allocator) :
        StochasticPropagator(input, ham, rank_allocator) {}
