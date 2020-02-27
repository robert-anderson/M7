//
// Created by rja on 27/02/2020.
//

#include "StochasticPropagator.h"

StochasticPropagator::StochasticPropagator(const std::unique_ptr <Hamiltonian> &ham,
                                           double tau, defs::ham_comp_t shift) : Propagator(ham, tau, shift) {}

