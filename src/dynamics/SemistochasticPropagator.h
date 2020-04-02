//
// Created by rja on 27/02/2020.
//

#ifndef M7_SEMISTOCHASTICPROPAGATOR_H
#define M7_SEMISTOCHASTICPROPAGATOR_H

#include <bits/unique_ptr.h>
#include "StochasticPropagator.h"
#include "DeterministicPropagator.h"

#if 0

class SemistochasticPropagator : public StochasticPropagator {
    std::unique_ptr<DeterministicPropagator> m_detprop = nullptr;
public:
    SemistochasticPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                             const RankAllocator<Determinant> &rank_allocator);

};


#endif //M7_SEMISTOCHASTICPROPAGATOR_H
#endif //M7_SEMISTOCHASTICPROPAGATOR_H
