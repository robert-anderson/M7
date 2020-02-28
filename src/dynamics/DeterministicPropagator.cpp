//
// Created by Robert John Anderson on 2020-02-11.
//

#include "DeterministicPropagator.h"

DeterministicPropagator::DeterministicPropagator(const std::unique_ptr<Hamiltonian> &ham,
                                                 const RankAllocator<Determinant> &rankAllocator,
                                                 defs::ham_comp_t target_norm, double tau,
                                                 defs::ham_comp_t shift) :
        Propagator(ham, rankAllocator, target_norm, tau, shift) {}
