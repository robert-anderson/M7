//
// Created by rja on 27/02/2020.
//

#include <src/defs.h>
#include "SemistochasticPropagator.h"

SemistochasticPropagator::SemistochasticPropagator(const std::unique_ptr<Hamiltonian> &ham, double tau,
                                                   const defs::ham_comp_t &shift) : StochasticPropagator(ham, tau,
                                                                                                         shift) {}
