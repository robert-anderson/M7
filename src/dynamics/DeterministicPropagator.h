//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_DETERMINISTICPROPAGATOR_H
#define M7_DETERMINISTICPROPAGATOR_H

#include "Propagator.h"

class DeterministicPropagator : public Propagator {
public:
    DeterministicPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                                                     const RankAllocator<Determinant> &rank_allocator);
};


#endif //M7_DETERMINISTICPROPAGATOR_H
