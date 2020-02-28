//
// Created by rja on 27/02/2020.
//

#include <src/defs.h>
#include "SemistochasticPropagator.h"

SemistochasticPropagator::SemistochasticPropagator(const std::unique_ptr<Hamiltonian> &ham,
                                                   const RankAllocator<Determinant> &rank_allocator,
                                                   defs::ham_comp_t target_shift, double tau,
                                                   defs::ham_comp_t shift, size_t seed) :
        StochasticPropagator(ham, rank_allocator, target_shift, tau, shift, seed) {}
