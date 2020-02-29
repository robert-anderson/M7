//
// Created by Robert John Anderson on 2020-02-11.
//

#include "DeterministicPropagator.h"

DeterministicPropagator::DeterministicPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                                                 const RankAllocator<Determinant> &rank_allocator) :
        Propagator(input, ham, rank_allocator) {}
