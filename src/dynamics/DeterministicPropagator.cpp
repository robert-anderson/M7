//
// Created by Robert John Anderson on 2020-02-11.
//

#include "DeterministicPropagator.h"

DeterministicPropagator::DeterministicPropagator(const std::unique_ptr<Hamiltonian> &ham,
                                                 const RankAllocator<Determinant> &rankAllocator, double tau,
                                                 defs::ham_comp_t shift) : Propagator(ham, rankAllocator, tau, shift) {}
