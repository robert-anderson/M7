//
// Created by rja on 27/02/2020.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H


#include <src/sample/PRNG.h>
#include "Propagator.h"

class StochasticPropagator : public Propagator {
    PRNG m_prng;
public:
    StochasticPropagator(const std::unique_ptr<Hamiltonian> &ham, const RankAllocator<Determinant> &rankAllocator,
                         double tau, defs::ham_comp_t shift, size_t seed);
};


#endif //M7_STOCHASTICPROPAGATOR_H
