//
// Created by rja on 27/02/2020.
//

#include "StochasticPropagator.h"


StochasticPropagator::StochasticPropagator(const std::unique_ptr<Hamiltonian> &ham,
                                           const RankAllocator<Determinant> &rank_allocator,
                                           defs::ham_comp_t target_shift, double tau,
                                           defs::ham_comp_t shift, size_t seed) :
        Propagator(ham, rank_allocator, target_shift, tau, shift), m_prng(seed) {}