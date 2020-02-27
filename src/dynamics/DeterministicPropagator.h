//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_DETERMINISTICPROPAGATOR_H
#define M7_DETERMINISTICPROPAGATOR_H

#include "Propagator.h"

class DeterministicPropagator : public Propagator {

    DeterministicPropagator(const std::unique_ptr<Hamiltonian> &ham,
               double tau, defs::ham_comp_t shift);
};


#endif //M7_DETERMINISTICPROPAGATOR_H
