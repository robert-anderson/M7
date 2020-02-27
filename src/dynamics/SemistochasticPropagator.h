//
// Created by rja on 27/02/2020.
//

#ifndef M7_SEMISTOCHASTICPROPAGATOR_H
#define M7_SEMISTOCHASTICPROPAGATOR_H

#include <bits/unique_ptr.h>
#include "StochasticPropagator.h"
#include "DeterministicPropagator.h"

class SemistochasticPropagator : public StochasticPropagator{
    std::unique_ptr<DeterministicPropagator> m_detprop = nullptr;
public:
    SemistochasticPropagator(const std::unique_ptr <Hamiltonian> &ham, double tau, const defs::ham_comp_t &shift);

};


#endif //M7_SEMISTOCHASTICPROPAGATOR_H
