//
// Created by rja on 27/02/2020.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H


class StochasticPropagator : public Propagator {
    PRNG m_prng;
public:
    StochasticPropagator(const std::unique_ptr<Hamiltonian> &ham, double tau, defs::ham_comp_t shift);
};


#endif //M7_STOCHASTICPROPAGATOR_H
